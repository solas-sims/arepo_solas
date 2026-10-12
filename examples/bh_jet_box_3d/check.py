""" @package ./examples/bh_jet_box_3d/check.py
Sanity checks for the BH jet box. Not a reference-solution comparison; it tests that the jet ran and that the books close.

usage: python3 check.py <simulation_directory> [True|False for plots]

Checks, using the snapshots (full precision) and log.txt (write the run's stdout there: mpirun ... > log.txt):
  1. total energy gain between snapshots equals the jet energy taken up by the gas, to 1e-3 relative (1e-2 with the beta-model
     potential, where the energy includes the potential energy of the gas, read from the constants in Config.sh)
  2. the change in net momentum is below 1% of the sum of |p| (the lobes are launched symmetrically, but the kick depends on each
     cell's own momentum, so it is not exactly zero)
  3. the jet tracer is set, and only around the launch radius / jet axis
"""

import re
import sys
import numpy as np
import h5py

simulation_directory = str(sys.argv[1])
makeplots = not (len(sys.argv) > 2 and sys.argv[2] != "True")

print("examples/bh_jet_box_3d/check.py: checking simulation output in directory " + simulation_directory)

LaunchRadius = 5.0  # kpc; BhJetLaunchRadius in param.txt

""" energy given/taken by the jet, per sync point, from the log """
t = None
taken = {}
for line in open(simulation_directory + "/log.txt"):
    m = re.match(r"Sync-Point \d+, Time: ([-+\d.eE]+)", line)
    if m:
        t = float(m.group(1))
    m = re.search(r"Energy given by BHs = ([-+\d.eE]+), Energy taken up by gas particles = ([-+\d.eE]+)", line)
    if m and t is not None:
        taken[t] = float(m.group(2))
times = sorted(taken)


def energy_taken_by(t_snap):
    prior = [x for x in times if x <= t_snap * (1 + 1e-9)]
    return taken[prior[-1]] if prior else 0.0


def read(isnap):
    f = h5py.File(simulation_directory + "/output/snap_%03d.hdf5" % isnap, "r")
    g = f["PartType0"]
    d = dict(
        time=f["Header"].attrs["Time"],
        pos=g["Coordinates"][:],
        vel=g["Velocities"][:],
        mass=g["Masses"][:],
        u=g["InternalEnergy"][:],
        box=f["Header"].attrs["BoxSize"],
    )
    if "PassiveScalars" in g:
        d["tracer"] = g["PassiveScalars"][:][:, 0] if g["PassiveScalars"].ndim > 1 else g["PassiveScalars"][:]
    f.close()
    return d


""" beta-model external potential (STATICBETA): Phi = 1.5 beta (k T / mu m_p) ln(1 + r^2/r_c^2), centred on the box centre """
beta_model = None
try:
    cfg = open(simulation_directory + "/Config.sh").read()
    if any(l.split("#")[0].strip() == "STATICBETA" for l in cfg.splitlines()):
        def macro(name):
            for l in cfg.splitlines():
                l = l.split("#")[0].strip()
                if l.startswith(name + "="):
                    return float(l.split("=")[1])
        UL, UV = 3.086e21, 1.0e5   # unit length [cm] and velocity [cm/s] of param.txt
        a = 1.38065e-16 * macro("BETA_MODEL_TEMP") / (macro("BETA_MODEL_MU") * 1.67262178e-24) / UV**2
        beta_model = (macro("BETA_MODEL_BETA"), macro("BETA_MODEL_RC"), a)
except OSError:
    pass


def etot(d):
    e = np.sum(d["mass"] * (d["u"] + 0.5 * np.sum(d["vel"] ** 2, axis=1)))
    if beta_model is not None:
        b, rc, a = beta_model
        c = np.array([0.5 * d["box"], 0.5 * d["box"], d["box"]])  # LONG_Z = 2
        r2 = np.sum((d["pos"] - c) ** 2, axis=1)
        e += np.sum(d["mass"] * 1.5 * b * a * np.log(1.0 + r2 / rc**2))
    return e


first = read(0)
isnap = 1
while True:
    try:
        last = read(isnap)
        isnap += 1
    except OSError:
        break
if isnap == 1:
    print("no snapshot after the IC found")
    sys.exit(1)

ok = True

dE = etot(last) - etot(first)
tk = energy_taken_by(last["time"])
rel = abs(dE - tk) / max(abs(tk), 1e-300)
# With the external potential the sum also carries the energy errors of the gravity coupling and of evaluating the potential at the
# mesh-generating points (a run without the jet drifts by a similar amount), so the tolerance is looser there.
tol = 1e-2 if beta_model is not None else 1e-3
print("t = %g: total energy gain %.8g, jet energy taken up %.8g (relative difference %.2e, tolerance %.0e)" % (last["time"], dE, tk, rel, tol))
if tk <= 0 or rel > tol:
    print("FAILED: energy does not close")
    ok = False

# change in net momentum since the IC (a relaxed IC has some net momentum already)
p = np.sum(last["mass"][:, None] * last["vel"], axis=0) - np.sum(first["mass"][:, None] * first["vel"], axis=0)
pabs = np.sum(last["mass"] * np.sqrt(np.sum(last["vel"] ** 2, axis=1)))
print("change in net momentum %s, sum |p| = %.4g" % (p, pabs))
if pabs > 0 and np.max(np.abs(p)) / pabs > 1e-2:
    print("FAILED: the net momentum changed by more than 1 per cent of sum |p|")
    ok = False

if "tracer" in last:
    jet = last["tracer"] > 0.01
    centre = np.array([0.5 * last["box"], 0.5 * last["box"], last["box"]])  # LONG_Z = 2
    r = np.sqrt(np.sum((last["pos"][jet] - centre) ** 2, axis=1))
    print("%d cells carry jet tracer" % np.sum(jet))
    if np.sum(jet) == 0:
        print("FAILED: no tracer")
        ok = False
    else:
        print("their distance from the BH: min %.2f, median %.2f, max %.2f kpc (launch radius %.1f)" % (r.min(), np.median(r), r.max(), LaunchRadius))

if makeplots:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        sel = np.abs(last["pos"][:, 0] - 0.5 * last["box"]) < 0.03 * last["box"]  # slab through the BH, about one cell thick
        fig, ax = plt.subplots(1, 2, figsize=(9, 6))
        speed = np.sqrt(np.sum(last["vel"] ** 2, axis=1))
        for a, q, name in zip(ax, [np.log10(speed + 1), last.get("tracer", speed)], ["log10(1+|v|) [km/s]", "jet tracer"]):
            s = a.scatter(last["pos"][sel, 1], last["pos"][sel, 2], c=q[sel], s=4)
            a.set_xlabel("y [kpc]")
            a.set_ylabel("z [kpc]")
            a.set_aspect("equal")
            plt.colorbar(s, ax=a, label=name, orientation="horizontal")
        fig.suptitle("t = %g" % last["time"])
        fig.savefig(simulation_directory + "/jet_slice.png", dpi=120)
        print("wrote jet_slice.png")
    except ImportError:
        print("matplotlib not available, skipping the plot")

sys.exit(0 if ok else 1)
