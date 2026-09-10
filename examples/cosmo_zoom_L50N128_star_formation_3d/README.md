# Cosmological box: 10 Mpc/h, 128^3, cooling + star formation + feedback (3D)

A small periodic cosmological box (10 Mpc/h, 128^3 dark matter particles,
Planck-ish cosmology) run from z=99 to z=4 (`TimeBegin=0.01`,
`TimeMax=0.2`) with `COOLING`+`USE_GRACKLE`, `AGORA_SF` star formation,
`WINDS`+`SUPERNOVAE` feedback, and on-the-fly `FOF`+`SUBFIND`. Gas is not
in the IC file -- it's split off from the dark matter at startup via
`GENERATE_GAS_IN_ICS`, based on `OmegaBaryon`.

Like `agora_disc_star_formation_3d`, this example does **not** plug into
`test.sh` -- it downloads external ICs and data tables rather than
generating them with `create.py`. Run it by hand, following the steps
below. `get_ics.sh` and `get_tables.sh` write to `./ICs/`,
`./grackle_data/`, and `./star_tables/` relative to the current working
directory, so run them from inside this directory (or wherever you've
copied it).

## What's in this directory

| File | Purpose |
|---|---|
| `Config.sh` | Compile-time flags |
| `param.txt` | Runtime parameters |
| `get_ics.sh` | Downloads the initial conditions into `./ICs/` |
| `get_tables.sh` | Downloads the Grackle cooling table and stellar feedback table into `./grackle_data/` and `./star_tables/` |
| `check.py` | Post-run sanity checks (see below -- not a reference-solution comparison) |

## 1. Get the initial conditions

```bash
./get_ics.sh
```

Fetches `cosmo_L10N128.hdf5` (~64 MB) from Dropbox into `./ICs/`. Not
committed to the repo (gitignored, same `examples/*/ICs/` rule as the
AGORA example).

The file is dark-matter-only: `PartType1` with 2,097,152 particles
(128^3), `BoxSize=10`, `Redshift=99` (`Time=0.01`), and cosmology
(`Omega0=0.3121`, `OmegaLambda=0.6879`, `HubbleParam=0.6751`) matching
`param.txt` exactly -- verified directly against the file's `/Header`.
`Config.sh` enables `GENERATE_GAS_IN_ICS`, which splits each DM particle
into a DM+gas pair at startup according to `OmegaBaryon`; the code
terminates on startup if the IC already contains gas, so don't feed it a
file that has `PartType0`.

Softening lengths (`0.003` Mpc/h, comoving) were carried over unchanged
from a previous L20/N256 run this param.txt was adapted from -- valid
here because both boxes have the same particle number density
(20/256 = 10/128 = 0.078125 Mpc/h mean interparticle spacing), so the
same softening-to-spacing ratio applies.

## 2. Get the data tables

```bash
./get_tables.sh
```

Downloads (gitignored, same pattern as `ICs/`):

- `./grackle_data/CloudyData_UVB=HM2012_high_density.h5` -- Grackle
  cooling/UV-background table, from
  [grackle_data_files](https://github.com/grackle-project/grackle_data_files)
  (pinned to commit `9286964`). Different table than the AGORA disc
  example uses (`FG2011_shielded.h5`) -- this one was specified for this
  run.
- `./star_tables/star_feedback_tables.hdf5` -- same stellar feedback
  lookup table as `agora_disc_star_formation_3d` (required by
  `WINDS`+`SUPERNOVAE` -> `STAR_FEEDBACK_ACTIVE`); see that example's
  README for the HDF5 schema `src/stars/star_tables.c` expects.

`TreecoolFile` (required by `COOLING`) is **not** downloaded -- it ships
with the repo at [`data/TREECOOL_ep`](../../data/TREECOOL_ep), which
`param.txt` already points at via `../../data/TREECOOL_ep`.

## 3. Build

You'll need a Grackle install (`GRACKLE_DIR` in `Makefile.systype`) in
addition to the usual MPI/GSL/HDF5/FFTW requirements described in
[`documentation/source/running.md`](../../documentation/source/running.md)
-- `PMGRID=256` in `Config.sh` means FFTW is required here, unlike the
AGORA disc example.

```bash
make CONFIG=examples/cosmo_L10N128_star_formation_3d/Config.sh \
     BUILD_DIR=examples/cosmo_L10N128_star_formation_3d/build \
     EXEC=examples/cosmo_L10N128_star_formation_3d/Arepo
```

## 4. Run

```bash
mkdir -p examples/cosmo_L10N128_star_formation_3d/output
cd examples/cosmo_L10N128_star_formation_3d
mpiexec -np <N> ./Arepo ./param.txt
```

## 5. Check

```bash
python check.py . True
```

No reference solution exists for this example yet -- nobody has run it
to completion to generate one. `check.py` is a smoke test rather than a
correctness check: per-snapshot it verifies the dark matter particle
count is unchanged from the IC, there are no NaN/Inf or negative masses,
positions stay inside the box, and (on the last snapshot) that star
particles actually formed; from `sfr.txt` it checks the cumulative
stellar mass is positive and monotonically non-decreasing, and that
total system mass (DM + gas + stars) hasn't drifted by more than 1%
between the first and last snapshot. Pass `True`/`False` as the second
argument to control whether it also writes plots to `./plots/`.

## Known gaps / things to revisit

- No reference-solution comparison, unlike `cosmo_box_star_formation_3d`
  -- once someone runs this to completion, consider saving a reduced
  `sfr.txt` as a committed reference and switching `check.py` to compare
  against it, the way the toy cosmo box example does.
- Build/run not verified locally (no Grackle install / MPI cluster in
  the environment this was written in).
- `param.txt`'s star formation and feedback parameters
  (`NumberDensThreshold=0.1 cm^-3`, `TemperatureThreshold=1e6 K`,
  `StarFormationEfficiency=0.5`, `SN_LeadTime=40 Myr`,
  `SN_HostShellSweepFrac=0.1`, `IMF=1` for Chabrier,
  `InitMetallicityinSolar=1e-10`) were supplied as-is and not
  independently re-derived or validated here.
- `check.py` doesn't inspect the FOF/SUBFIND group catalogs
  (`groups_*/fof_subhalo_tab_*`) even though `Config.sh` enables them --
  only particle-level snapshot data and `sfr.txt` are checked.
