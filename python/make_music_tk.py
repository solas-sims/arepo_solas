#!/usr/bin/env python3
"""Write a MUSIC `transfer = music` transfer-function table from CAMB or CLASS.

Output columns (whitespace separated, NO header, as read by transfer_music.cc):
    k[h/Mpc]  T_tot  T_cdm  T_baryon  T_v_cdm  T_v_baryon

Conventions chosen here (MUSIC takes log10 of every column, wants T -> const at
low k, and replaces values <= 0 by a k^-2 extrapolation):
  * k in h/Mpc
  * densities divided by k^2 where needed, signs flipped so low-k is positive
  * all columns divided by one common factor so T_tot(kmin) = 1 (ratios kept)
    -> MUSIC normalises with sigma_8 itself, so the overall amplitude is free.

Velocity columns:  --velocities placeholder  (default) copies the cdm/baryon
density columns. That is only safe for DM-only MUSIC runs, where MUSIC does not
use them. --velocities real writes genuine velocity transfers, whose
normalisation convention versus MUSIC is NOT verified (see README text printed
at the end); check MUSIC's input_powerspec_vcdm.txt before trusting baryon ICs.

Parameters may come from a YAML file (--params cosmo.yaml); keys are the option
names (code, h, omega_m, omega_b, ns, zstart, As, mnu, neff, kmax, velocities,
output, k_per_decade, nk). Command-line options override the file. Needs PyYAML.

Example:
  python make_music_tk.py camb  --h 0.6774 --omega_b 0.0486 --omega_m 0.3089 \
      --ns 0.9667 --zstart 49 -o tk_camb.dat
  python make_music_tk.py class --h 0.6774 --omega_b 0.0486 --omega_m 0.3089 \
      --ns 0.9667 --zstart 49 -o tk_class.dat
"""
import argparse
import sys

import numpy as np


# ----------------------------------------------------------------------------
def run_camb(a):
    import camb
    from camb import model

    h = a.h
    ombh2 = a.omega_b * h**2
    omch2 = (a.omega_m - a.omega_b) * h**2
    pars = camb.set_params(H0=100 * h, ombh2=ombh2, omch2=omch2, ns=a.ns,
                           As=a.As, mnu=a.mnu, num_massive_neutrinos=1 if a.mnu > 0 else 0,
                           nnu=a.neff, omk=0.0, tau=0.06)
    kw = {"k_per_logint": max(1, int(round(a.k_per_decade / np.log(10.0))))} if a.k_per_decade else {}
    pars.set_matter_power(redshifts=[a.zstart], kmax=a.kmax, **kw)
    pars.NonLinear = model.NonLinear_none
    res = camb.get_results(pars)
    td = res.get_matter_transfer_data()
    data = td.transfer_data  # shape (nvars, nk, nz); index = model.Transfer_* - 1

    def col(name):
        return np.array(data[getattr(model, name) - 1, :, 0], dtype=float)

    k = col("Transfer_kh")  # h/Mpc
    out = dict(k=k, tot=col("Transfer_tot"), cdm=col("Transfer_cdm"), b=col("Transfer_b"))
    if a.velocities == "real":
        out["vcdm"] = col("Transfer_Newt_vel_cdm")
        out["vb"] = col("Transfer_Newt_vel_baryon")
    return out


# ----------------------------------------------------------------------------
def run_class(a):
    from classy import Class

    h = a.h
    base = {
        "h": h,
        "omega_b": a.omega_b * h**2,
        "omega_cdm": (a.omega_m - a.omega_b) * h**2,
        "n_s": a.ns,
        "A_s": a.As,
        "N_ur": (a.neff - 1.0132) if a.mnu > 0 else a.neff,  # one massive species ~1.0132 of Neff
        "z_pk": a.zstart,
        "P_k_max_h/Mpc": a.kmax,
        "output": "mTk,vTk",
    }
    if a.k_per_decade:
        base["k_per_decade_for_pk"] = a.k_per_decade
    if a.mnu > 0:
        base["N_ncdm"] = 1
        base["m_ncdm"] = a.mnu

    def get(gauge):
        c = Class()
        c.set(dict(base, gauge=gauge))
        c.compute()
        tr = c.get_transfer(z=a.zstart)
        c.struct_cleanup()
        c.empty()
        return tr

    # densities: synchronous gauge (same convention as CAMB's CDM/baryon transfers)
    s = get("synchronous")
    k = np.array(s["k (h/Mpc)"], dtype=float)
    k2 = k**2
    # CLASS d_x ~ -T k^2 at low k: divide by k^2 to get a plateau; sign fixed later
    out = dict(k=k, tot=np.array(s["d_tot"]) / k2, cdm=np.array(s["d_cdm"]) / k2,
               b=np.array(s["d_b"]) / k2)
    if a.velocities == "real":
        # theta_cdm vanishes identically in synchronous gauge -> use Newtonian gauge
        n = get("newtonian")
        kn = np.array(n["k (h/Mpc)"], dtype=float)
        if kn.shape != k.shape or not np.allclose(kn, k):
            raise RuntimeError("CLASS k grids differ between gauges; cannot combine")
        out["vcdm"] = np.array(n["t_cdm"]) / k2
        out["vb"] = np.array(n["t_b"]) / k2
    return out


# ----------------------------------------------------------------------------
def resample(d, nk):
    """Cubic spline in log k / log T onto nk log-spaced points (columns are positive here;
    non-positive columns are splined linearly in T instead). Adds no new physics."""
    from scipy.interpolate import CubicSpline
    lk = np.log10(d["k"])
    new_lk = np.linspace(lk[0], lk[-1], nk)
    out = {"k": 10.0**new_lk}
    for key, v in d.items():
        if key == "k":
            continue
        if np.all(v > 0):
            out[key] = 10.0 ** CubicSpline(lk, np.log10(v))(new_lk)
        else:
            out[key] = CubicSpline(lk, v)(new_lk)
    return out


def finalise(d, a):
    k = d["k"]
    order = np.argsort(k)
    d = {key: v[order] for key, v in d.items()}
    keys = ["tot", "cdm", "b"] + (["vcdm", "vb"] if a.velocities == "real" else [])

    # make each column positive at lowest k
    for key in keys:
        if d[key][0] < 0:
            d[key] = -d[key]
    # common normalisation: T_tot(kmin) = 1 (keeps relative amplitudes)
    f = d["tot"][0]
    for key in keys:
        d[key] = d[key] / f
    if a.velocities != "real":
        d["vcdm"], d["vb"] = d["cdm"].copy(), d["b"].copy()

    if np.any(np.diff(d["k"]) <= 0):
        raise RuntimeError("k is not strictly increasing")
    if a.nk:
        d = resample(d, a.nk)
    for key in ["tot", "cdm", "b", "vcdm", "vb"]:
        bad = np.count_nonzero(d[key] <= 0)
        if bad:
            print(f"WARNING: column '{key}' has {bad}/{len(d['k'])} non-positive values; "
                  f"MUSIC will replace them by a k^-2 extrapolation (baryon "
                  f"oscillations / sign conventions are the usual cause).",
                  file=sys.stderr)
    return d


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("code", nargs="?", choices=["camb", "class"], help="or set `code:` in the YAML file")
    p.add_argument("--params", help="YAML file of parameters (keys = option names below); "
                                    "command-line options override it")
    p.add_argument("--h", type=float)
    p.add_argument("--omega_m", type=float, help="total matter (cdm+b), not times h^2")
    p.add_argument("--omega_b", type=float, help="baryons, not times h^2")
    p.add_argument("--ns", type=float)
    p.add_argument("--zstart", type=float, help="must equal MUSIC [setup] zstart")
    p.add_argument("--As", type=float, default=2.1e-9, help="irrelevant for MUSIC (sigma_8 renormalises)")
    p.add_argument("--mnu", type=float, default=0.0, help="sum m_nu in eV (one massive species)")
    p.add_argument("--neff", type=float, default=3.046)
    p.add_argument("--kmax", type=float, default=500.0, help="h/Mpc")
    p.add_argument("--k_per_decade", type=float, default=None,
                   help="ask CAMB (k_per_logint) / CLASS (k_per_decade_for_pk) for finer native k sampling")
    p.add_argument("--nk", type=int, default=None,
                   help="resample the final table onto this many log-spaced k points (needs SciPy)")
    p.add_argument("--velocities", choices=["placeholder", "real"], default="placeholder")
    p.add_argument("-o", "--output")

    # YAML values become defaults; anything given on the command line wins.
    pre, _ = p.parse_known_args()
    if pre.params:
        import yaml
        with open(pre.params) as fh:
            y = yaml.safe_load(fh) or {}
        known = {act.dest for act in p._actions} - {"help", "params"}
        unknown = set(y) - known
        if unknown:
            p.error(f"unknown key(s) in {pre.params}: {sorted(unknown)}; valid: {sorted(known)}")
        p.set_defaults(**y)
    a = p.parse_args()
    missing = [n for n in ("code", "h", "omega_m", "omega_b", "ns", "zstart", "output")
               if getattr(a, n) is None]
    if missing:
        p.error(f"missing required value(s): {missing} (give on command line or in YAML)")
    if a.code not in ("camb", "class"):
        p.error("code must be 'camb' or 'class'")
    if a.velocities not in ("placeholder", "real"):
        p.error("velocities must be 'placeholder' or 'real'")

    d = run_camb(a) if a.code == "camb" else run_class(a)
    d = finalise(d, a)

    tab = np.column_stack([d["k"], d["tot"], d["cdm"], d["b"], d["vcdm"], d["vb"]])
    np.savetxt(a.output, tab, fmt="%.8e")  # no header: MUSIC would parse '#' lines as data
    print(f"wrote {len(tab)} rows to {a.output}  (k = {d['k'][0]:.3e} .. {d['k'][-1]:.3e} h/Mpc, z = {a.zstart})")
    print("MUSIC config: [cosmology] transfer = music ; transfer_file = <this file> ; "
          "set zstart, sigma_8 there.")
    if a.velocities == "placeholder":
        print("NOTE: velocity columns are placeholders (copies of density) -> DM-only runs only.")
    else:
        print("NOTE: real velocity columns: convention vs MUSIC unverified; compare MUSIC's "
              "input_powerspec_vcdm.txt / _vbaryon.txt before using for baryon ICs.")


if __name__ == "__main__":
    main()
