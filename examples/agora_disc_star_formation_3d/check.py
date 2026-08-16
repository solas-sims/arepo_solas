""" @package ./examples/agora_disc_star_formation_3d/check.py
Sanity checks and monitoring across the snapshots of an
agora_disc_star_formation_3d run.

Not a reference-solution comparison (no golden run exists for this
example -- it depends on downloaded external data, see get_ics.sh /
get_tables.sh) and not wired into test.sh for the same reason. Run by
hand after a run:

    python check.py <simulation_directory>

where <simulation_directory> contains output/snap_*.hdf5.

Two pass/fail invariant checks:

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

Plus three monitoring checks, printed per-snapshot (mostly informational,
but each still fails on an obviously-broken run):

  Star formation
    Number of star particles per snapshot. Stars are never removed from
    the particle list once formed in this code (only "deactivated" when
    dead/low-mass -- see deactivate_star() in src/stars/star_update.c),
    so this must be non-decreasing.

  Star-forming gas
    Mass and mass fraction of gas above the AGORA_SF density threshold
    and below the temperature threshold (src/star_formation/sfr_AGORA.c's
    sf_criteria()), i.e. eligible to form stars. Temperature is computed
    with the same non-Grackle mean-molecular-weight fallback formula
    sf_criteria() itself uses when USE_GRACKLE is off (mu = 4/(8-5*(1-XH)),
    fully-ionized approximation) -- this is an approximation of the
    actual (Grackle-chemistry-network-dependent) trigger condition used
    at runtime, not an exact reproduction; treat as directional/trend
    monitoring, not a precise reconstruction of the code's own decision.
    Thresholds are read from each snapshot's own embedded /Parameters
    group, not hardcoded.

  Stellar spatial extent
    Min/median/max cylindrical radius and |z|-height of star particles
    relative to the box centre. Flagged (not hard-failed) if the maximum
    radius approaches the box edge, which would indicate particles have
    been numerically ejected rather than staying in the disc.

Exits 0 if both pass/fail checks pass (dust check may be skipped), 1
otherwise. The three monitoring checks print WARN rather than FAIL for
anything short of a clear invariant violation, since "trend looks odd"
isn't itself proof of a bug the way a metals/dust invariant violation is.
"""

import glob
import os
import sys

import h5py
import numpy as np

FloatType = np.float64

PROTONMASS_CGS = 1.67262178e-24
BOLTZMANN_CGS = 1.38065e-16
HYDROGEN_MASSFRAC = 0.76
GAMMA_MINUS1 = 2.0 / 3.0  # monatomic gas, gamma = 5/3


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


def check_star_formation(snapshot_files):
    print("--- star formation check ---")
    counts = []
    times = []
    for filename in snapshot_files:
        with h5py.File(filename, "r") as f:
            n_stars = int(f["Header"].attrs["NumPart_Total"][4])
            counts.append(n_stars)
            times.append(float(f["Header"].attrs["Time"]))
            print("  t=%.6g  number of stars=%d" % (times[-1], n_stars))

    ok = True
    for i in range(1, len(counts)):
        if counts[i] < counts[i - 1]:
            print(
                "FAIL: number of stars decreased from %d (t=%.6g) to %d (t=%.6g) "
                "-- stars should never be removed from the particle list"
                % (counts[i - 1], times[i - 1], counts[i], times[i])
            )
            ok = False

    if ok:
        print("PASS: number of stars is non-decreasing across all snapshots")
    return ok


def sf_thresholds(f):
    """Read NumberDensThreshold/TemperatureThreshold from the snapshot's own
    embedded /Parameters group (not hardcoded). Returns None if not present
    (e.g. the run wasn't built with AGORA_SF)."""
    if "Parameters" not in f:
        return None
    attrs = f["Parameters"].attrs
    if "NumberDensThreshold" not in attrs or "TemperatureThreshold" not in attrs:
        return None
    return float(attrs["NumberDensThreshold"]), float(attrs["TemperatureThreshold"])


def gas_number_density_and_temperature(f):
    """Approximate physical number density [cm^-3] and temperature [K] per
    gas cell, using the same formulas as sf_criteria() in
    src/star_formation/sfr_AGORA.c, but with the non-Grackle mean molecular
    weight fallback (mu = 4/(8-5*(1-XH))) regardless of whether the run
    actually used USE_GRACKLE -- an approximation of the true
    chemistry-network-dependent mu, not an exact match. See module
    docstring."""
    header = f["Header"].attrs
    unit_length_cm = float(header["UnitLength_in_cm"])
    unit_mass_g = float(header["UnitMass_in_g"])
    unit_velocity_cm_s = float(header["UnitVelocity_in_cm_per_s"])
    unit_density_cgs = unit_mass_g / unit_length_cm ** 3

    density = np.array(f["PartType0"]["Density"], dtype=FloatType)
    utherm = np.array(f["PartType0"]["InternalEnergy"], dtype=FloatType)

    mu = 4.0 / (8.0 - 5.0 * (1.0 - HYDROGEN_MASSFRAC))

    number_density = (density * unit_density_cgs) / mu / PROTONMASS_CGS
    temperature = (
        utherm * unit_velocity_cm_s ** 2 * mu * PROTONMASS_CGS * GAMMA_MINUS1 / BOLTZMANN_CGS
    )

    return number_density, temperature


def check_sf_gas_fraction(snapshot_files):
    print("--- star-forming gas check (approximate, see module docstring) ---")
    with h5py.File(snapshot_files[0], "r") as f:
        thresholds = sf_thresholds(f)
        if thresholds is None:
            print("SKIP: no NumberDensThreshold/TemperatureThreshold in /Parameters "
                  "(run was not built with AGORA_SF)")
            return True
    number_dens_threshold, temperature_threshold = thresholds
    print(
        "  thresholds from snapshot: NumberDensThreshold=%.6g cm^-3, TemperatureThreshold=%.6g K"
        % (number_dens_threshold, temperature_threshold)
    )

    ok = True
    for filename in snapshot_files:
        with h5py.File(filename, "r") as f:
            time = float(f["Header"].attrs["Time"])
            mass = np.array(f["PartType0"]["Masses"], dtype=FloatType)
            number_density, temperature = gas_number_density_and_temperature(f)

            eligible = (number_density >= number_dens_threshold) & (temperature <= temperature_threshold)
            eligible_mass = float(np.sum(mass[eligible]))
            total_mass = float(np.sum(mass))
            fraction = eligible_mass / total_mass if total_mass > 0 else float("nan")

            if not np.isfinite(fraction):
                print("FAIL: t=%.6g: star-forming gas fraction is not finite" % time)
                ok = False

            print(
                "  t=%.6g  star-forming-eligible gas mass=%.6e (%.4g%% of gas mass), cells=%d/%d"
                % (time, eligible_mass, 100.0 * fraction, int(np.sum(eligible)), len(mass))
            )

    if ok:
        print("PASS: star-forming gas fraction computed successfully at every snapshot "
              "(informational -- no expected trend asserted)")
    return ok


def check_stellar_extent(snapshot_files):
    print("--- stellar spatial extent check ---")
    ok = True
    for filename in snapshot_files:
        with h5py.File(filename, "r") as f:
            time = float(f["Header"].attrs["Time"])
            boxsize = float(np.asarray(f["Header"].attrs["BoxSize"]).reshape(-1)[0])
            centre = boxsize / 2.0

            if "PartType4" not in f or int(f["Header"].attrs["NumPart_Total"][4]) == 0:
                print("  t=%.6g  no star particles yet" % time)
                continue

            pos = np.array(f["PartType4"]["Coordinates"], dtype=FloatType)
            dxy = pos[:, :2] - centre
            r = np.hypot(dxy[:, 0], dxy[:, 1])
            z = np.abs(pos[:, 2] - centre)

            print(
                "  t=%.6g  R [kpc]: min=%.3g median=%.3g max=%.3g | |z| [kpc]: min=%.3g median=%.3g max=%.3g"
                % (time, r.min(), np.median(r), r.max(), z.min(), np.median(z), z.max())
            )

            # loose sanity bound: flag (not hard-fail) if stars approach the
            # box edge, which would indicate numerical ejection rather than
            # staying in the disc
            edge_fraction = 0.9
            if r.max() > edge_fraction * (boxsize / 2.0) or z.max() > edge_fraction * (boxsize / 2.0):
                print(
                    "WARN: t=%.6g: a star particle is within %.0f%% of the box edge "
                    "(R_max=%.3g, |z|_max=%.3g, box half-size=%.3g) -- check for numerical ejection"
                    % (time, 100.0 * edge_fraction, r.max(), z.max(), boxsize / 2.0)
                )

    print("PASS: stellar spatial extent computed successfully at every snapshot with stars "
          "(informational, WARN only on apparent ejection)")
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
    stars_ok = check_star_formation(snapshot_files)
    sf_gas_ok = check_sf_gas_fraction(snapshot_files)
    extent_ok = check_stellar_extent(snapshot_files)

    if metals_ok and dust_ok and stars_ok and sf_gas_ok and extent_ok:
        sys.exit(0)
    else:
        sys.exit(1)
