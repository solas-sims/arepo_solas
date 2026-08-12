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

## What's in this directory

| File | Purpose |
|---|---|
| `Config.sh` | Compile-time flags (copy of the root `Config_AGORA.sh`, moved here) |
| `param.txt` | Runtime parameters, derived from the values baked into the IC file itself (see below) |
| `get_ics.sh` | Downloads the initial conditions into `./ICs/` |
| `get_tables.sh` | Downloads the Grackle cooling table and stellar feedback table into `./grackle_data/` and `./star_tables/` |
| `agora_disc_diagnostics.py` | Post-processing: radial surface density / scale height / vertical velocity dispersion profiles, plus a smoothed density-grid projection, from one or more output snapshots |

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

### `./grackle_data/CloudyData_UVB=HM2012_high_density.h5`

Grackle cooling/UV-background table, from
[grackle_data_files](https://github.com/grackle-project/grackle_data_files)
(pinned to commit `9286964`). Matches `GRACKLE_CHEMISTRY=3` in `Config.sh`.

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
    --output disc_profiles.png --grid-output "disc_grid_{name}.png"
```

Produces overlaid radial surface-density / scale-height / vertical velocity
dispersion profiles across the given snapshots, plus one smoothed
density-grid projection image per snapshot. Uses the group's `analysistools`
package if installed, otherwise falls back to a self-contained `h5py`-based
reader — see the script's docstring for the full set of options
(`--mask-radius`, `--ngrid`, `--save-data`, etc).

## Known gaps / things to revisit

- `param.txt`'s star formation and feedback parameters
  (`NumberDensThreshold=5 cm^-3`, `TemperatureThreshold=1e4 K`,
  `StarFormationEfficiency=0.01`, `SN_LeadTime=3 Myr`,
  `SN_HostShellSweepFrac=0.5`, `IMF=0` for Kroupa) are literature-typical
  defaults, not values validated against this fork's implementation —
  treat as a starting point.
- No `check.py`/reference data, unlike the `test.sh`-integrated examples —
  there's nothing here yet to catch a silent regression automatically.
