""" @package ./examples/agora_disc_star_formation_3d/check.py
Sanity checks on metal and dust conservation/invariants across the
snapshots of an agora_disc_star_formation_3d run.

Not a reference-solution comparison (no golden run exists for this
example -- it depends on downloaded external data, see get_ics.sh /
get_tables.sh) and not wired into test.sh for the same reason. Run by
hand after a run:

    python check.py <simulation_directory>

where <simulation_directory> contains output/snap_*.hdf5.

Two checks:

  Metals check
    Total metal mass (gas-phase + whatever's locked into star particles)
    must never decrease over the run, within a small floating-point
    tolerance -- metals are created by stellar enrichment and never
    destroyed anywhere in this code.

  Dust check (only if the run was built with DUST; skipped with a
  warning otherwise)
    Every gas cell's DustMass must satisfy the module's core physical
    invariant at every snapshot: 0 <= DustMass <= (gas-phase metal mass
    in that cell). Dust is a subset of the metal budget it condensed
    from, so it can never be negative or exceed it.

Exits 0 if both checks pass (or the dust check is skipped because the
run has no DustMass field), 1 otherwise.
"""

import glob
import os
import sys

import h5py
import numpy as np

FloatType = np.float64


def find_snapshots(simulation_directory):
    directory = os.path.join(simulation_directory, "output")
    files = sorted(glob.glob(os.path.join(directory, "snap_*.hdf5")))
    if not files:
        print("check.py: no snap_*.hdf5 files found in " + directory)
        sys.exit(1)
    return files


def gas_metal_mass(f):
    """Metal mass per gas cell (METALS_INDEX == 0 in PassiveScalars)."""
    mass = np.array(f["PartType0"]["Masses"], dtype=FloatType)
    passive_scalars = np.array(f["PartType0"]["PassiveScalars"], dtype=FloatType)
    metallicity = passive_scalars[:, 0]
    return metallicity * mass, mass


def star_metal_mass(f):
    """Metal mass locked into star particles, if any exist yet."""
    if "PartType4" not in f or "Metallicity" not in f["PartType4"]:
        return 0.0
    mass = np.array(f["PartType4"]["Masses"], dtype=FloatType)
    metallicity = np.array(f["PartType4"]["Metallicity"], dtype=FloatType)
    return float(np.sum(metallicity * mass))


def check_metals(snapshot_files):
    print("--- metals check ---")
    total_metal_mass = []
    times = []
    for filename in snapshot_files:
        with h5py.File(filename, "r") as f:
            metal_mass_per_cell, _ = gas_metal_mass(f)
            total = float(np.sum(metal_mass_per_cell)) + star_metal_mass(f)
            total_metal_mass.append(total)
            times.append(float(f["Header"].attrs["Time"]))

    total_metal_mass = np.array(total_metal_mass)
    times = np.array(times)

    ok = True
    # allow a tiny relative tolerance for floating-point noise, but total
    # metal mass must not meaningfully decrease between snapshots
    rel_tol = 1e-6
    for i in range(1, len(total_metal_mass)):
        prev, cur = total_metal_mass[i - 1], total_metal_mass[i]
        if cur < prev * (1.0 - rel_tol):
            print(
                "FAIL: total metal mass decreased from %.6e (t=%.6g) to %.6e (t=%.6g)"
                % (prev, times[i - 1], cur, times[i])
            )
            ok = False

    for t, m in zip(times, total_metal_mass):
        print("  t=%.6g  total metal mass=%.6e" % (t, m))

    if ok:
        print("PASS: total metal mass is non-decreasing across all snapshots")
    return ok


def check_dust(snapshot_files):
    print("--- dust check ---")
    with h5py.File(snapshot_files[0], "r") as f:
        if "DustMass" not in f["PartType0"]:
            print("SKIP: no DustMass field in snapshots (run was not built with DUST)")
            return True

    ok = True
    tol = 1e-20  # absolute tolerance, code mass units
    for filename in snapshot_files:
        with h5py.File(filename, "r") as f:
            time = float(f["Header"].attrs["Time"])
            dust_mass = np.array(f["PartType0"]["DustMass"], dtype=FloatType)
            metal_mass_per_cell, _ = gas_metal_mass(f)

            n_negative = int(np.sum(dust_mass < -tol))
            n_exceeds = int(np.sum(dust_mass > metal_mass_per_cell + tol))

            if n_negative > 0:
                print("FAIL: t=%.6g: %d cells have DustMass < 0" % (time, n_negative))
                ok = False
            if n_exceeds > 0:
                print(
                    "FAIL: t=%.6g: %d cells have DustMass > GasMetals (dust exceeds its own metal budget)"
                    % (time, n_exceeds)
                )
                ok = False

            print(
                "  t=%.6g  total dust mass=%.6e  nonzero cells=%d/%d"
                % (time, float(np.sum(dust_mass)), int(np.sum(dust_mass > 0)), len(dust_mass))
            )

    if ok:
        print("PASS: 0 <= DustMass <= GasMetals holds in every cell, every snapshot")
    return ok


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: python check.py <simulation_directory>")
        sys.exit(1)

    simulation_directory = str(sys.argv[1])
    print(
        "examples/agora_disc_star_formation_3d/check.py: checking simulation output in directory "
        + simulation_directory
    )

    snapshot_files = find_snapshots(simulation_directory)

    metals_ok = check_metals(snapshot_files)
    dust_ok = check_dust(snapshot_files)

    if metals_ok and dust_ok:
        sys.exit(0)
    else:
        sys.exit(1)
