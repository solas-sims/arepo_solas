# AGORA isolated disc + star formation (3D)

An isolated Milky-Way-mass disc galaxy (gas + stellar disc + bulge + dark
matter halo), run with the AGORA-style star formation model (`AGORA_SF`),
Grackle cooling/chemistry, and this fork's stellar feedback module
(`STARS`/`STAR_PARTICLES`/`WINDS`/`SUPERNOVAE`). This closes the gap noted in
[`documentation/source/solas-parameters.md`](../../documentation/source/solas-parameters.md)
that no committed example turns on `AGORA_SF`.

This example does **not** plug into `test.sh` — it downloads external ICs
and data tables rather than generating them with `create.py`. Run it by
hand, following the steps below.

`get_ics.sh` and `get_tables.sh` write to `./ICs/`, `./grackle_data/`, and
`./star_tables/` relative to the current working directory (not to wherever
the scripts themselves live), so run them from inside this directory — or
from wherever you've copied it, e.g. a `test.sh`-style run directory.

## What's in this directory

| File | Purpose |
|---|---|
| `Config.sh` | Compile-time flags (copy of the root `Config_AGORA.sh`, moved here) |
| `param.txt` | Runtime parameters, derived from the values baked into the IC file itself (see below) |
| `get_ics.sh` | Downloads the initial conditions into `./ICs/` |
| `get_tables.sh` | Downloads the Grackle cooling table and stellar feedback table into `./grackle_data/` and `./star_tables/` |
| `agora_disc_diagnostics.py` | Post-processing: radial surface density / scale height / vertical velocity dispersion profiles, plus a smoothed density-grid projection, from one or more output snapshots |
| `check.py` | Sanity checks (not a reference-solution comparison): total metal mass never decreases across snapshots, and every gas cell's `DustMass` stays within `0 <= DustMass <= GasMetals` at every snapshot (skipped with a message if the run wasn't built with `DUST`) |

## 1. Get the initial conditions

```bash
./get_ics.sh
```

This fetches `agora_lowres_arepo_ic_without_u-with-grid.hdf5` (~37 MB) from
Dropbox into `./ICs/`. It is **not** committed to the repo (see
`.gitignore`'s `examples/*/ICs/` rule) — re-run the script instead of
committing the file.

This particular file is the AGORA low-resolution disc IC *after* the
`ADDBACKGROUNDGRID` step (mesh already added at `GridSize=32`), not the raw
particle IC — it's ready to feed straight into an Arepo physics run. It was
itself produced by someone else's Arepo run (visible in its embedded
`/Config` and `/Parameters` HDF5 groups — `Git_date` "Wed Jan 21 11:24:30
2026 +0800", output path under `/scratch/pawsey1164/nmaragkakis/GALAXY/`).
`param.txt` in this directory copies the physical setup (box size,
units, softenings, mesh regularization settings) straight out of those
embedded groups so the two stay consistent:

- `BoxSize = 1200` (kpc), non-cosmological (`ComovingIntegrationOn=0`,
  `PeriodicBoundariesOn=0`)
- Units: `UnitLength_in_cm=3.08568e21` (1 kpc), `UnitMass_in_g=1.989e42`
  (1e9 M_sun — note this is **not** the same mass unit as
  `examples/galaxy_merger_star_formation_3d`, which uses 1e10 M_sun),
  `UnitVelocity_in_cm_per_s=1e5` (1 km/s)
- `NumPart_Total = [752516, 100000, 100000, 12500, 0, 0]` — gas, DM halo,
  stellar disc, stellar bulge
- `ReferenceGasPartMass = 8.59322e-05` (in code units, ≈ 8.6e4 M_sun)
- All gas cells have `InternalEnergy = 0` in the file (hence the
  `without_u` in the filename) — `InitGasTemp` in `param.txt` sets the
  actual starting temperature (10000 K, adjust as needed)

## 2. Get the data tables

```bash
./get_tables.sh
```

Downloads two files (gitignored, same pattern as `ICs/`):

### `./grackle_data/CloudyData_UVB=FG2011_shielded.h5`

Grackle cooling/UV-background table, from
[grackle_data_files](https://github.com/grackle-project/grackle_data_files)
(pinned to commit `9286964`).

### `./star_tables/star_feedback_tables.hdf5`

Required because `Config.sh` enables `WINDS` + `SUPERNOVAE`, which the
Makefile combines into `STAR_FEEDBACK_ACTIVE`
([`Makefile:442`](../../Makefile)) — a stellar-evolution lookup table
(mass-loss rates, wind velocities, SN ejecta, as a function of stellar mass
and metallicity), read by
[`src/stars/star_tables.c`](../../src/stars/star_tables.c). Verified against
the schema reverse-engineered from `star_tables.c`: `Z_COUNT=13`
(`Z_VALUES` from `1e-11` to `0.03`), `M_COUNT=32` (`M_VALUES` from `2` to
`120` M_sun), and every `Z=.../M=.../` group carries `Age`, `Radius`,
`Temperature`, `MassLossRate`, `HLossRate`, `HeLossRate`, `MetalsLossRate`,
`WindVelocity`, `SN_MassLoss`, `SN_HLoss`, `SN_HeLoss`, `SN_MetalsLoss` (plus
`Energy`/`Photons`, unused here since `RADIATION` is off).

Units for these fields aren't documented anywhere in `src/stars/` — worth
double-checking with whoever built the table generation pipeline if
feedback strengths look physically off.

### `TreecoolFile` (COOLING)

Unlike the two tables above, this one ships with the repo already — no
download needed. `param.txt`'s `TreecoolFile` points at
`../../TREECOOL_ep`, i.e. the repo root's copy; an identical copy also
lives at [`data/TREECOOL_ep`](../../data/TREECOOL_ep) in the Arepo source
directory if you'd rather reference that one.

## 3. Build

You'll need a Grackle install (`GRACKLE_DIR` in `Makefile.systype`, see
`Makefile`'s `GRACKLE_DIR ?= $(HOME)/software/grackle_solas/` default) in
addition to the usual MPI/GSL/HDF5 requirements described in
[`documentation/source/running.md`](../../documentation/source/running.md).

```bash
make CONFIG=examples/agora_disc_star_formation_3d/Config.sh \
     BUILD_DIR=examples/agora_disc_star_formation_3d/build \
     EXEC=examples/agora_disc_star_formation_3d/Arepo
```

## 4. Run

```bash
mkdir -p examples/agora_disc_star_formation_3d/output
cd examples/agora_disc_star_formation_3d
mpiexec -np <N> ./Arepo ./param.txt
```

## 5. Diagnostics

```bash
python agora_disc_diagnostics.py --snapshot output/snap_000.hdf5 output/snap_010.hdf5 \
    --output disc_profiles.png --grid-output "disc_grid_{name}.png" \
    --metals-grid-output "disc_grid_metals_{name}.png" \
    --dust-grid-output "disc_grid_dust_{name}.png" \
    --dtg-grid-output "disc_grid_dust_to_gas_{name}.png" \
    --metals-dust-output "disc_metals_dust_profiles.png"
```

Produces overlaid radial surface-density / scale-height / vertical velocity
dispersion profiles across the given snapshots, plus one smoothed
density-grid projection image per snapshot. Uses the group's `analysistools`
package if installed, otherwise falls back to a self-contained `h5py`-based
reader — see the script's docstring for the full set of options
(`--mask-radius`, `--ngrid`, `--save-data`, etc).

Grid projections use each gas cell's own smoothing length (effective radius
derived from `Density`/`Masses`, i.e. an SPH-kernel-style deposition) rather
than a single fixed smoothing scale for the whole box — the AGORA IC's cell
size varies by orders of magnitude across the disc, so one global scale either
over-smooths the dense inner region or leaves particle-scale gaps in the
sparse outer disc. Falls back to a fixed `--smooth-sigma` (CIC + Gaussian)
only when a snapshot has no `Density` field, or `--ptype` isn't gas. The
colourbar is a proper column density (Σ, summed along the projection axis and
divided by pixel area — not a mean of raw per-cell values), converted to
M☉ pc⁻² using each snapshot's own `UnitMass_in_g`. All three panels (XY/XZ/YZ)
of a given plot share one colour scale by default (computed from that plot's
own combined data range); pass `--grid-vmin`/`--grid-vmax` to fix it
explicitly, e.g. for comparing the same panel across multiple snapshots.

Also reads each snapshot's gas-cell metal mass (`PassiveScalars`) and, if
present, dust mass (`DustMass`, only written when `Config.sh` has `DUST`
enabled) directly from `PartType0`. This produces the same smoothed-grid
XY/XZ/YZ projection as gas mass gets (`disc_grid_metals_*`/
`disc_grid_dust_*`), plus a dust-to-gas mass ratio grid
(`disc_grid_dust_to_gas_*` — gas and dust mass deposited and projected
independently, then divided; dimensionless, so unaffected by the M☉ pc⁻²
conversion) and a second overlaid-profile figure
(`disc_metals_dust_profiles.png`) of metal and dust surface density vs.
radius, mirroring the gas surface-density panel. The dust/dust-to-gas
grids and profile are skipped (with a log message, not an error) if the
run has no `DustMass` field.

## 6. Sanity checks

```bash
python check.py .
```

Run from this directory (or pass whatever directory contains your
`output/`). Five checks, printed per-snapshot and summarized pass/fail:

- **Metals**: total metal mass (gas-phase + whatever's locked into star
  particles) must never decrease across the run — metals are only ever
  created by stellar enrichment in this code, never destroyed.
- **Dust** (only meaningful if `Config.sh` has `DUST` enabled; skipped
  with a message otherwise): every gas cell's `DustMass` must satisfy
  `0 <= DustMass <= GasMetals` at every snapshot — dust is a subset of
  the metal budget it condensed from.
- **Star formation**: number of star particles must never decrease —
  stars are only "deactivated" when dead, never removed from the
  particle list.
- **Star-forming gas** (informational, no hard pass/fail beyond "did it
  compute"): mass and mass fraction of gas above the `AGORA_SF` density
  threshold and below the temperature threshold, using thresholds read
  from each snapshot's own embedded `/Parameters` group. Temperature
  uses the same non-Grackle mean-molecular-weight fallback the C code
  itself uses when `USE_GRACKLE` is off — an approximation of the real
  (Grackle-chemistry-dependent) trigger condition, not an exact
  reproduction, so treat this as trend monitoring rather than a precise
  reconstruction of what actually decided each cell's eligibility.
- **Stellar spatial extent** (informational, `WARN` rather than `FAIL`):
  min/median/max cylindrical radius and `|z|`-height of star particles
  relative to the box centre, flagged if a star approaches the box edge
  (a sign of numerical ejection rather than staying in the disc).

Not a reference-solution comparison like the `test.sh`-integrated
examples' `check.py` scripts — there's no golden run for this example
(external IC/data dependencies), so these are invariant/sanity checks
only, meant to catch an obviously broken run (mass appearing from
nowhere, dust exceeding its metal budget) rather than verify exact
physics.

## Known gaps / things to revisit

- `InitMetallicityinSolar = 0.1` (required by `METALS`) is an assumed
  starting value, not one taken from the IC file or validated against this
  fork's implementation.
- `param.txt`'s star formation and feedback parameters
  (`NumberDensThreshold=5 cm^-3`, `TemperatureThreshold=1e4 K`,
  `StarFormationEfficiency=0.01`, `SN_LeadTime=3 Myr`,
  `SN_HostShellSweepFrac=0.5`, `IMF=0` for Kroupa) are literature-typical
  defaults, not values validated against this fork's implementation —
  treat as a starting point.
