""" @package ./examples/cosmo_L10N128_star_formation_3d/check.py
Sanity checks for the L10N128 cosmological box run (cooling, AGORA_SF,
WINDS+SUPERNOVAE feedback, FOF/SUBFIND).

Unlike the other cosmo_box_*_3d examples, this one starts from a
downloaded, externally-generated IC (see get_ics.sh) rather than one
built by create.py, and no simulation has been run against it yet to
produce reference output. So this script does not compare against a
reference solution -- it checks that the run did something physically
sane: mass conservation, no NaN/Inf/negative masses, dark matter particle
count unchanged, star formation actually happened, gas temperatures in a
sane range. Treat it as a smoke test, not a correctness proof.
"""

import sys
import os
import glob
import numpy as np
import h5py

simulation_directory = str(sys.argv[1])
print("cosmo_L10N128_star_formation_3d: checking simulation output in directory " + simulation_directory)

makeplots = True
if len(sys.argv) > 2:
    makeplots = sys.argv[2] == "True"

FloatType = np.float64

directory = simulation_directory + "/output/"

# ---- reference quantities from the IC / param.txt ----
BoxSize = 10.0  # Mpc/h
HubbleParam = 0.6751
NumPartDM_IC = 128 ** 3  # 2097152, dark-matter-only IC; gas is split off via GENERATE_GAS_IN_ICS

errors = []


def fail(msg):
    errors.append(msg)
    print("FAIL: " + msg)


def find_snapshots():
    files = sorted(glob.glob(directory + "snap_*.hdf5"))
    if not files:
        fail("no snap_*.hdf5 files found in " + directory)
    return files


def check_snapshot(path, is_first, is_last, initial_total_mass):
    with h5py.File(path, "r") as f:
        header = f["Header"].attrs
        boxsize = float(np.asarray(header["BoxSize"]).reshape(-1)[0])
        if not np.isclose(boxsize, BoxSize, rtol=1e-6):
            fail(f"{path}: BoxSize={boxsize}, expected {BoxSize}")

        num_part = np.asarray(header["NumPart_Total"], dtype=np.int64)

        # dark matter count (PartType1) must never change -- no accretion/merging in this Config
        if num_part[1] != NumPartDM_IC:
            fail(f"{path}: NumPart_Total[1] (DM) = {num_part[1]}, expected {NumPartDM_IC}")

        total_mass = 0.0
        for ptype in range(6):
            key = f"PartType{ptype}"
            if key not in f or num_part[ptype] == 0:
                continue
            grp = f[key]

            if "Masses" in grp:
                mass = np.asarray(grp["Masses"], dtype=FloatType)
            else:
                mass_table = np.asarray(header["MassTable"], dtype=FloatType)
                mass = np.full(num_part[ptype], mass_table[ptype], dtype=FloatType)

            if np.any(~np.isfinite(mass)):
                fail(f"{path}: PartType{ptype} has non-finite Masses")
            if np.any(mass < 0):
                fail(f"{path}: PartType{ptype} has negative Masses (min={mass.min()})")
            total_mass += float(np.sum(mass))

            if "Coordinates" in grp:
                pos = np.asarray(grp["Coordinates"], dtype=FloatType)
                if np.any(~np.isfinite(pos)):
                    fail(f"{path}: PartType{ptype} has non-finite Coordinates")
                if np.any(pos < -1e-6) or np.any(pos > boxsize + 1e-6):
                    fail(f"{path}: PartType{ptype} has Coordinates outside [0, BoxSize]")

            if "Velocities" in grp:
                vel = np.asarray(grp["Velocities"], dtype=FloatType)
                if np.any(~np.isfinite(vel)):
                    fail(f"{path}: PartType{ptype} has non-finite Velocities")

            if ptype == 0 and "InternalEnergy" in grp:
                u = np.asarray(grp["InternalEnergy"], dtype=FloatType)
                if np.any(~np.isfinite(u)):
                    fail(f"{path}: PartType0 has non-finite InternalEnergy")
                if np.any(u < 0):
                    fail(f"{path}: PartType0 has negative InternalEnergy")

        if is_last:
            if num_part[4] == 0:
                fail(f"{path}: no star particles (PartType4) formed by the end of the run -- "
                     "check NumberDensThreshold/TemperatureThreshold/StarFormationEfficiency in param.txt")

        return total_mass, num_part


def check_sfr_log():
    path = directory + "sfr.txt"
    if not os.path.exists(path):
        fail(path + " not found (expected from USE_SFR)")
        return

    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    cum_mass_stars = data[:, 5]
    if np.any(~np.isfinite(cum_mass_stars)):
        fail("sfr.txt: non-finite cumulative stellar mass")
    if np.any(np.diff(cum_mass_stars) < -1e-8):
        fail("sfr.txt: cumulative stellar mass is not monotonically non-decreasing")
    if cum_mass_stars[-1] <= 0:
        fail("sfr.txt: cumulative stellar mass is zero at the end of the run -- no star formation occurred")

    sfr_rate = data[:, 3]
    if np.any(sfr_rate < -1e-8):
        fail("sfr.txt: negative star formation rate")

    return data


sfr_data = check_sfr_log()

snapshots = find_snapshots()
initial_total_mass = None
final_total_mass = None
if snapshots:
    for i, path in enumerate(snapshots):
        total_mass, num_part = check_snapshot(
            path, is_first=(i == 0), is_last=(i == len(snapshots) - 1), initial_total_mass=initial_total_mass
        )
        if i == 0:
            initial_total_mass = total_mass
        if i == len(snapshots) - 1:
            final_total_mass = total_mass

    if initial_total_mass is not None and final_total_mass is not None:
        rel_change = abs(final_total_mass - initial_total_mass) / initial_total_mass
        # winds/SNe recycle mass between gas and stars, so this is a loose bound,
        # not a tight conservation check
        if rel_change > 0.01:
            fail(f"total mass changed by {rel_change:.2%} between first and last snapshot "
                 f"(initial={initial_total_mass:g}, final={final_total_mass:g})")

if makeplots and snapshots:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plots_dir = simulation_directory + "/plots"
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)

    if sfr_data is not None:
        fig = plt.figure()
        ax = plt.axes([0.15, 0.15, 0.8, 0.8])
        ax.plot(sfr_data[:, 0], sfr_data[:, 5])
        ax.set_xlabel("scale factor")
        ax.set_ylabel(r"cumulative stellar mass [code units]")
        ax.set_yscale("log")
        fig.savefig(plots_dir + "/cumulative_stellar_mass.pdf")
        plt.close(fig)

    with h5py.File(snapshots[-1], "r") as f:
        pos_dm = np.asarray(f["PartType1"]["Coordinates"], dtype=FloatType)
        fig = plt.figure(figsize=(6, 6))
        ax = plt.axes([0.12, 0.12, 0.85, 0.85])
        if pos_dm.shape[0] > 50000:
            sel = np.random.choice(pos_dm.shape[0], 50000, replace=False)
        else:
            sel = np.arange(pos_dm.shape[0])
        ax.scatter(pos_dm[sel, 0], pos_dm[sel, 1], marker=".", s=0.1, alpha=0.3, rasterized=True)
        header = f["Header"].attrs
        num_part = np.asarray(header["NumPart_Total"], dtype=np.int64)
        if num_part[4] > 0:
            pos_stars = np.asarray(f["PartType4"]["Coordinates"], dtype=FloatType)
            ax.scatter(pos_stars[:, 0], pos_stars[:, 1], marker="*", c="r", s=1.0, alpha=0.8, rasterized=True)
        ax.set_xlim([0, BoxSize])
        ax.set_ylim([0, BoxSize])
        ax.set_xlabel("[Mpc/h]")
        ax.set_ylabel("[Mpc/h]")
        fig.savefig(plots_dir + "/dark_matter_and_stars_final.pdf", dpi=200)
        plt.close(fig)

if errors:
    print(f"\n{len(errors)} check(s) failed.")
    sys.exit(1)

print("all checks passed.")
sys.exit(0)
