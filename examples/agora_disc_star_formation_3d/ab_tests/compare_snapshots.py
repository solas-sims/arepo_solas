#!/usr/bin/env python3
"""Compare two Arepo output directories snapshot by snapshot (A/B test of batched vs per-cell Grackle).

usage: compare_snapshots.py DIR_A DIR_B [--rtol 0.0] [--snaps 0 1 2]

Particles are matched by ParticleIDs. For every dataset in every PartType group the script reports whether the values are
bit-identical and otherwise the largest absolute and relative differences. Exit status is 1 if any difference exceeds rtol
(default 0, i.e. require bit-identical output) or if particle sets differ.
"""
import argparse, glob, os, re, sys
import h5py
import numpy as np


def snaps(d):
    out = {}
    for f in glob.glob(os.path.join(d, "snap_*.hdf5")) + glob.glob(os.path.join(d, "snapdir_*", "snap_*.0.hdf5")):
        m = re.search(r"snap_(\d+)", os.path.basename(f))
        if m:
            out[int(m.group(1))] = f
    return out


def load(path):
    with h5py.File(path, "r") as f:
        data = {}
        for g in f:
            if not g.startswith("PartType"):
                continue
            ids = f[g]["ParticleIDs"][...]
            order = np.argsort(ids, kind="stable")
            data[g] = {k: f[g][k][...][order] for k in f[g] if isinstance(f[g][k], h5py.Dataset)}
        t = f["Header"].attrs["Time"]
    return t, data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--rtol", type=float, default=0.0)
    ap.add_argument("--snaps", type=int, nargs="*")
    args = ap.parse_args()

    sa, sb = snaps(args.a), snaps(args.b)
    common = sorted(set(sa) & set(sb))
    if args.snaps is not None:
        common = [s for s in common if s in args.snaps]
    if not common:
        sys.exit("no common snapshots found")
    if set(sa) ^ set(sb):
        print("WARNING: snapshots only in one run:", sorted(set(sa) ^ set(sb)))

    bad = False
    for s in common:
        ta, da = load(sa[s])
        tb, db = load(sb[s])
        print(f"== snapshot {s:03d}  time A={ta:.8g} B={tb:.8g}")
        for g in sorted(set(da) | set(db)):
            if g not in da or g not in db:
                print(f"  {g}: present in only one run"); bad = True; continue
            na, nb = len(da[g]["ParticleIDs"]), len(db[g]["ParticleIDs"])
            if na != nb or not np.array_equal(da[g]["ParticleIDs"], db[g]["ParticleIDs"]):
                print(f"  {g}: particle sets differ (N_A={na}, N_B={nb})"); bad = True; continue
            for k in sorted(da[g]):
                if k not in db[g]:
                    continue
                x, y = da[g][k], db[g][k]
                if np.array_equal(x, y, equal_nan=True) if x.dtype.kind == "f" else np.array_equal(x, y):
                    continue
                if x.dtype.kind not in "fiu":
                    continue
                d = np.abs(x.astype(float) - y.astype(float))
                rel = d / np.maximum(np.abs(x.astype(float)), np.abs(y.astype(float)), 1e-300)
                nd = int(np.count_nonzero(d))
                flag = "" if rel.max() <= args.rtol else "  <-- exceeds rtol"
                bad |= bool(flag)
                print(f"  {g}/{k}: {nd}/{x.size} differ, max abs {d.max():.3e}, max rel {rel.max():.3e}{flag}")
        print("  (datasets not listed are bit-identical)")
    print("RESULT:", "DIFFERENCES FOUND" if bad else "OK")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
