# Cosmological zoom: Halo 381 in a 50 Mpc/h parent box, cooling + star formation + feedback (3D)

A zoom-in on one halo (dubbed "Halo 381") selected from a periodic L50/h
Mpc, N128 parent box (Planck-ish cosmology), run from z=99 to z=4
(`TimeBegin=0.01`, `TimeMax=0.2`) with `COOLING`+`USE_GRACKLE`,
`AGORA_SF` star formation, `WINDS`+`SUPERNOVAE` feedback, and on-the-fly
`FOF`+`SUBFIND`. The high-res region is refined to `levelmax=9` (an
effective 512^3 grid); everything outside it is the coarse `levelmin=7`
(128^3) boundary population. Gas is not in the IC file -- it's split off
from the dark matter at startup via `GENERATE_GAS_IN_ICS`, based on
`OmegaBaryon`.

Like `agora_disc_star_formation_3d` and `cosmo_L10N128_star_formation_3d`,
this example does **not** plug into `test.sh` -- it generates ICs with an
external tool (MUSIC2) and downloads external data tables rather than
using `create.py`. Run it by hand, following the steps below. `get_ics.sh`
writes to `./ICs/` and `get_tables.sh` writes to `./grackle_data/` and
`./star_tables/`, both relative to the current working directory, so run
them from inside this directory (or wherever you've copied it).

## How Halo 381 was selected

Zoom ICs need a Lagrangian patch (`ref_center`/`ref_extent`) identified
in advance, by running a uniform parent box, picking a halo from its
z=0 FOF catalogue, and tracing that halo's particles back to the parent
ICs. That's already been done for this example -- `get_ics.sh` (step 1
below) does not redo it, it only runs MUSIC2 on the already-selected
patch. The provenance, for reproducing or picking a different halo:

1. `MUSIC2_conf/zoom_parent.conf` generates the unigrid L50/h Mpc, N128
   parent-box ICs (`levelmin=levelmax=7`).
2. That parent box was run through Arepo (gravity-only, to z=0) with FOF
   enabled, producing a snapshot + halo catalogue.
3. `MUSIC2_scripts/get_zoom_ics_from_parent.py` reads that snapshot and
   catalogue, selects a halo by mass range, traces its particles back to
   the parent ICs (via `analysistools`), and reports the Lagrangian
   patch centre/extent -- which is what's hardcoded as `ref_center`/
   `ref_extent` in `MUSIC2_conf/zoom_halo381.conf`.

To target a different halo, rerun step 3 against your own parent-box
run/catalogue with a different `--mass-min`/`--mass-max`/`--halo-index`,
and copy the reported centre/extent into a new MUSIC2 conf (set
`levelmax=9`, `levelmin_TF=9`, matching `zoom_halo381.conf`).

## What's in this directory

| File | Purpose |
|---|---|
| `Config.sh` | Compile-time flags |
| `param.txt` | Runtime parameters |
| `MUSIC2_conf/zoom_parent.conf` | MUSIC2 config for the unigrid parent box (halo-selection provenance) |
| `MUSIC2_conf/zoom_halo381.conf` | MUSIC2 config for the actual zoom ICs used by this run |
| `MUSIC2_scripts/get_zoom_ics_from_parent.py` | Halo selection / Lagrangian-patch tracing script |
| `get_ics.sh` | Runs MUSIC2 on `zoom_halo381.conf`, writing the zoom ICs into `./ICs/` |
| `get_tables.sh` | Downloads the Grackle cooling table and stellar feedback table into `./grackle_data/` and `./star_tables/` |
| `check.py` | Post-run sanity checks (see below -- not a reference-solution comparison) |

## 1. Get the initial conditions

```bash
./get_ics.sh
```

Runs MUSIC2 against `MUSIC2_conf/zoom_halo381.conf`, writing
`./ICs/L50_N128_z99_Halo381.hdf5` (dark-matter-only: `PartType1`
high-res + `PartType2` coarse boundary, `BoxSize=50`, `Redshift=99`
(`Time=0.01`), cosmology `Omega0=0.3121`, `OmegaLambda=0.6879`,
`HubbleParam=0.6751` matching `param.txt`). Needs a MUSIC2 build --
by default the script looks for one at `../../../MUSIC2/build/MUSIC`
(i.e. a sibling `MUSIC2` checkout next to this repo); point `MUSIC2_EXE`
at your build if it's elsewhere. Not committed to the repo (gitignored,
same `examples/*/ICs/` rule as the AGORA example).

`Config.sh` enables `GENERATE_GAS_IN_ICS`, which splits each DM particle
into a DM+gas pair at startup according to `OmegaBaryon`; the code
terminates on startup if the IC already contains gas, so don't feed it a
file that has `PartType0`.

Softening lengths in `param.txt` are set per particle type from the
level-7/level-9 grid spacing using MUSIC's own `boxsize/2^level/40`
convention (the same formula it prints as `suggested_highressoft` in the
IC header): `0.00244140625` Mpc/h for the level-9 high-res population
(gas/DM/stars, types 0/1/4) and `0.009765625` Mpc/h for the level-7
coarse boundary DM (type 2, matching `Config.sh`'s
`INDIVIDUAL_GRAVITY_SOFTENING=4`). These are formulaic, not
independently tuned or validated against the actual clustering in this
patch -- see "Known gaps" below.

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
-- `PMGRID=256` plus `PLACEHIGHRESREGION`/`PM_ZOOM_OPTIMIZED` in
`Config.sh` mean FFTW is required here.

```bash
make CONFIG=examples/cosmo_zoom_L50N128_star_formation_3d/Config.sh \
     BUILD_DIR=examples/cosmo_zoom_L50N128_star_formation_3d/build \
     EXEC=examples/cosmo_zoom_L50N128_star_formation_3d/Arepo
```

## 4. Run

```bash
mkdir -p examples/cosmo_zoom_L50N128_star_formation_3d/output
cd examples/cosmo_zoom_L50N128_star_formation_3d
mpiexec -np <N> ./Arepo ./param.txt
```

## 5. Check

```bash
python check.py . True
```

No reference solution exists for this example yet -- nobody has run it
to completion to generate one. `check.py` is a smoke test rather than a
correctness check: per-snapshot it verifies the high-res (`PartType1`)
and coarse boundary (`PartType2`) dark matter particle counts are
unchanged from the first snapshot, there are no NaN/Inf or negative
masses, positions stay inside the box, and (on the last snapshot) that
star particles actually formed; from `sfr.txt` it checks the cumulative
stellar mass is positive and monotonically non-decreasing, and that
total system mass (DM + gas + stars) hasn't drifted by more than 1%
between the first and last snapshot. Pass `True`/`False` as the second
argument to control whether it also writes plots to `./plots/`.

## Lessons carried over from the L10N128 box

These come from timing and consistency runs of
`cosmo_L10N128_star_formation_3d` on Setonix and are built into this
example's `Config.sh`/`param.txt`:

- **Star-particle timestep cap (`STAR_MAX_TIMESTEP_MYR=0.1`).** The
  hard-coded 0.01 Myr cap was written for resolved galaxy runs. In the
  box it drove the global step count to 177k (41.5 ks wall) and formed
  ~2.5x more stars than caps of 0.1 or 1 Myr, which agreed with each
  other to ~20% in stellar mass (0.1 Myr: 28.8k steps, 8.4 ks). 0.1 Myr
  is therefore used here. Single realisation per setting, so this is
  evidence of convergence, not proof. The option exists on the
  `grackle-batched-cooling-v2` branch (commit `e3b3e3c`); on a branch
  without it the line is inert and the old 0.01 Myr cap applies.
- **`STAR_FEEDBACK_IMMEDIATE_DEPOSIT`.** Changed the 0.1 Myr box result by
  under 1% in stellar mass (about 9% in SFR density) at no run-time
  cost, and removes any dependence on the non-conservative host-removal
  limiter, so it is on.
- **`GRACKLE_CHEMISTRY=1`.** The atomic network is enough for halos
  resolved here at z=4. H2/HD cooling (levels 2-3) matters in minihalos
  below ~1e8 Msun, which are not resolved. Snapshots and restart files
  are not compatible with chemistry-3 runs.
- **`REFINEMENT_HIGH_RES_GAS`** is specific to the zoom. It restricts
  refinement and derefinement to gas that came from the high-res region
  (`AllowRefinement`), carries a per-cell `HighResMass`, and gives a star
  formed in a cell with little high-res mass the coarse Type 3 softening.
  Here all gas is split from the level-9 `PartType1` particles at
  startup, so every cell starts as high-res and the flag changes little
  at the start. It adds `HighResGasMass` and `AllowRefinement` to the snapshots.
- **`DEBUG` and `EVALPOTENTIAL` off.** Neither is needed (no black hole
  seeding), and both cost run time.
- **Star-formation parameters.** The box example's committed
  `param.txt` said `StarFormationEfficiency=0.5`, `TemperatureThreshold=1e6 K`;
  the runs actually used 0.01 and 2e4 K. This example uses the latter.
  Even so the box overproduced stars at z=4 by roughly 10x relative to the
  Behroozi+2019 SMHM relation, so treat the star formation scheme
  (AGORA, fixed n and T) and feedback strength as unvalidated for
  cosmological runs.
- **Cooling batching.** Batched Grackle cooling gave a 2.3x faster
  `sfrcool` and bit-identical output to the per-cell path, but cooling is
  under 1% of total run time, so it does not shorten the run.
- **Build and run settings.** Build with `-O2`/`-O3`, not the Setonix
  `Makefile`'s default `-O0`. Raise `CpuTimeBetRestartFile` for
  production: restart writes were ~9% of the box's run time at 300 s.
  `param.txt` here still has 300, left unchanged.

## Known gaps / things to revisit

- No reference-solution comparison, unlike `cosmo_box_star_formation_3d`
  -- once someone runs this to completion, consider saving a reduced
  `sfr.txt` as a committed reference and switching `check.py` to compare
  against it, the way the toy cosmo box example does.
- Build/run not verified locally (no Grackle install / MPI cluster in
  the environment this was written in).
- Softening lengths (see step 1 above) follow MUSIC's `1/40 of grid
  spacing` convention per level but weren't independently re-derived or
  checked against this patch's actual particle spacing/clustering.
- Star formation and feedback parameters in `param.txt`
  (`StarFormationEfficiency=0.01`, `TemperatureThreshold=2e4 K`,
  `NumberDensThreshold=0.1 cm^-3`, `SN_LeadTime=40 Myr`,
  `SN_HostShellSweepFrac=0.1`, `IMF=1` for Chabrier,
  `InitMetallicityinSolar=1e-10`) are the values the 10 Mpc/h box runs
  used (see "Lessons carried over from the L10N128 box"). They are not
  tuned for this zoom: its gas cells are about twice as heavy
  (~1.9e7 vs ~0.95e7 Msun), and the density threshold in particular
  is not resolution-scaled. A scan of `NumberDensThreshold` at 1 and
  10 cm^-3 on the box was still running when this was written; revisit
  this value once it finishes.
- `check.py` doesn't inspect the FOF/SUBFIND group catalogs
  (`groups_*/fof_subhalo_tab_*`) even though `Config.sh` enables them --
  only particle-level snapshot data and `sfr.txt` are checked.
- The parent-box run (step in "How Halo 381 was selected") itself isn't
  scripted here -- only its MUSIC2 conf and the post-hoc analysis script
  are. Reproducing the halo selection from scratch means running that
  conf, building/running Arepo on it separately, and running
  `get_zoom_ics_from_parent.py` by hand.
