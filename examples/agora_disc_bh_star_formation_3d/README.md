# AGORA isolated disc + star formation + central black hole (3D)

An isolated Milky-Way-mass disc galaxy (gas + stellar disc + bulge + dark matter halo) with a
single central supermassive black hole. It is the
[`agora_disc_star_formation_3d`](../agora_disc_star_formation_3d/README.md) setup (AGORA-style star
formation `AGORA_SF`, Grackle cooling/chemistry, and this fork's stellar feedback module
`STARS`/`STAR_PARTICLES`/`WINDS`/`SUPERNOVAE`) plus the black-hole module:

- **BH in the ICs, no on-the-fly seeding.** A 1e6 Msun `PartType5` particle sits at the box
  centre (the AGORA discs are built about `BoxSize/2`), at rest in the construction frame.
  `BLACKHOLE_SEEDING`, `HALO_SEEDING` and FOF are *not* enabled.
- **Bondi accretion** (`BONDI_ACCRETION`), capped at the Eddington rate.
- **Thermal feedback only** (`BH_THERMAL_FEEDBACK`). Do not use `BH_FEEDBACK`: it also turns on the
  jet, which is not operational.
- Two resolutions, `param_lowres.txt` and `param_medres.txt`.

This example does **not** plug into `test.sh`: it downloads external ICs and data tables rather
than generating them with `create.py`. Run it by hand, following the steps below.

`get_ics.sh` and `get_tables.sh` write to `./ICs/`, `./grackle_data/` and `./star_tables/` relative
to the current working directory (not to wherever the scripts live), so run them from inside this
directory, or from wherever you have copied it.

## What's in this directory

| File | Purpose |
|---|---|
| `Config.sh` | Compile-time flags: the parent's, plus `BLACKHOLES`, `BONDI_ACCRETION`, `BH_THERMAL_FEEDBACK` |
| `param_lowres.txt`, `param_medres.txt` | Runtime parameters for the two resolutions (replace the parent's single `param.txt`) |
| `get_ics.sh` | Downloads the low- or medium-resolution IC into `./ICs/` and adds the central BH |
| `add_central_bh.py` | Appends one `PartType5` particle at `BoxSize/2` to an IC file (called by `get_ics.sh`) |
| `get_tables.sh` | Downloads the Grackle cooling table and the stellar feedback table (see section 2 for a caveat) |
| `agora_disc_diagnostics.py` | Post-processing: radial profiles and smoothed density-grid projections from output snapshots |

## What differs from the parent example

| | Parent | This example |
|---|---|---|
| `Config.sh` | -- | adds `BLACKHOLES`, `BONDI_ACCRETION`, `BH_THERMAL_FEEDBACK` |
| Parameters | -- | adds `BhDesNgb`, `BhDesDev`, `HMaxFactor`, `Epsilon_r`, `Epsilon_f`, `Mload` |
| IC | `agora_lowres_arepo_ic-with-grid.hdf5` | same file plus one `PartType5` (`add_central_bh.py`) |
| `param.txt` | one file | `param_lowres.txt`, `param_medres.txt` |

### Black hole parameters

| Parameter | Value | Meaning |
|---|---|---|
| `Epsilon_f` | 0.01 | fraction of the radiated luminosity coupled to the gas as thermal energy |
| `HMaxFactor` | 12.5 | maximum BH search radius = `HMaxFactor` x the BH's softening (type 5 -> softening type 1): 1.0 kpc low-res, 0.5 kpc med-res |
| `Epsilon_r` | 0.01 | radiative efficiency; sets the Eddington cap and the luminosity |
| `Mload` | 0 | mass-loading fraction. Energy deposited is `Epsilon_f (1 - Mload) Epsilon_r Mdot c^2`, mass returned `Mload Epsilon_r Mdot`, so 0 gives pure thermal feedback |
| `BhDesNgb`, `BhDesDev` | 32, 8 | smoothing-length targets, copied from the cosmological examples; untuned for a disc |

### Settings differing from the parent's `param.txt`

These come from the reference parameter files (`low.txt`, `med.txt`) supplied with this example,
not from the parent:

| Parameter | Parent | Here | Note |
|---|---|---|---|
| `InitMetallicityinSolar` | 0.1 | 1.0 | |
| `NumberDensThreshold` | 5 | 100 | cm^-3, physical density threshold |
| `IMF` | 0 (Kroupa) | 1 (Chabrier) | 2 = Salpeter |
| `SN_LeadTime` | 3 | 40 | Myr |
| `SN_HostShellSweepFrac` | 0.5 | 0.1 | internally capped at 0.9 |
| `GrackleDataFile` | `./grackle_data/CloudyData_UVB=FG2011_shielded.h5` | `./grackle_data/cloudy_metals_2008_3D.h5` | fetched by `get_tables.sh` |
| `TimeBetSnapshot` | 0.02 | 0.01 | ~9.8 Myr |
| `MaxSizeTimestep` | 0.0005 | 1e-2 | |
| `TimeBetStatistics` | 0.0234375 | 0 | |
| `TimeLimitCPU`, `CpuTimeBetRestartFile` | 93000, 12000 | 100000, 10000 | |
| `UnitLength_in_cm` | 3.08568e21 | 3.085678e21 | |
| `MeanVolume` | 52734.4 | 536377 | both resolutions |

`GrackleSpeciesState`, `RaySplitFactor` and `RayMaxDistance` are present in the param files for
consistency with the reference files but are not used by this build (`RADIATION` is off;
`GrackleSpeciesState` and `RayMaxDistance` are not read by this checkout at all).
`BoxSize` stays at each IC's header value (the reference files say 1300 for both).

### Resolution-dependent values

| | low-res | med-res |
|---|---|---|
| gas / halo / disc / bulge particles | 1e5 / 1e5 / 1e5 / 1.25e4 | 1e6 / 1e6 / 1e6 / 1.25e5 |
| gas particle mass | 8.6e4 Msun | 8.6e3 Msun |
| `BoxSize` | 1200 kpc | 1523.4424 kpc (from the IC) |
| softenings (all types) | 0.08 kpc | 0.04 kpc |
| `MinimumComovingHydroSoftening` | 0.08 | 0.04 |
| `ReferenceGasPartMass` | 8.59322e-05 | 8.59322e-06 |
| `MinVolume` | 5.12e-4 (dx >= 80 pc) | 6.40e-5 (dx >= 40 pc) |

Particle counts and gas masses are those of the parent's low-res IC and the previous medium-res IC;
they have not been re-checked against the current downloads.

Common to both (from the parent): non-cosmological (`ComovingIntegrationOn=0`,
`PeriodicBoundariesOn=0`); units 1 kpc, 1e9 Msun (`UnitMass_in_g=1.989e42`; not the 1e10 Msun of
`examples/galaxy_merger_star_formation_3d`) and 1 km/s; `TimeMax = 0.5` (~489 Myr).

## 1. Get the initial conditions

Requires `python3` with `h5py` and `numpy`. The ICs are **not** committed (`.gitignore`'s
`examples/*/ICs/` rule): re-run the script instead of committing the files.

### Low resolution

```bash
./get_ics.sh lowres
```

Downloads `agora_lowres_arepo_ic-with-grid.hdf5` from Dropbox into `./ICs/` and writes
`./ICs/agora_lowres_arepo_ic_with_bh.hdf5`. This is the AGORA low-resolution disc IC *after* the
`ADDBACKGROUNDGRID` step (mesh already added), so it is ready to feed into an Arepo physics run.
The gas `InternalEnergy` (~101 in code units) and the header `UnitMass_in_g` (1.989e42) are used
as supplied; the script only appends the BH.

### Medium resolution

```bash
./get_ics.sh medres [PATH_TO_IC]                 # writes ICs/agora_medres_arepo_ic_with_bh.hdf5
```

Downloads `agora_medres_arepo_ic-with-grid.hdf5` (the AGORA medium-resolution IC, also already
through `ADDBACKGROUNDGRID`) into `./ICs/`. To use a local copy instead, pass its path. As for low
resolution, units and gas `InternalEnergy` are used as supplied and the script only adds the BH,
at `BoxSize/2` = 761.72 kpc, the disc centre in this file. Other properties of the file:

| | medium-res file | note |
|---|---|---|
| disc / halo particle types | PartType1 = stellar disc, PartType2 = halo | **unchanged** (low-res is the reverse) |
| `BoxSize` | 1523.44 (embedded parameters say 2000) | `param_medres.txt` uses 1523.4424209594727 |

The types are not swapped in the file. Softenings are equal for types 1 to 3, so the code is not
affected, but anything assuming PartType2 is the disc (for example `agora_disc_diagnostics.py`)
must be pointed at PartType1.

`BH_MASS_MSUN=3e6 ./get_ics.sh ...` changes the BH mass (default 1e6 Msun). `add_central_bh.py`
gives the BH the next unused particle ID, updates the header particle counts, and prints the
offset of the stellar-disc centre of mass from the BH as a sanity check.

## 2. Get the data tables

```bash
./get_tables.sh
```

Downloads two files (gitignored, same pattern as `ICs/`):

### `./grackle_data/cloudy_metals_2008_3D.h5`

Grackle cooling table (Cloudy metals, 2008, 3D), downloaded from Dropbox. The param files point
`GrackleDataFile` at it. The parent example's `CloudyData_UVB=FG2011_shielded.h5` is no longer
downloaded here.

### `./star_tables/star_feedback_tables.hdf5`

Required because `Config.sh` enables `WINDS` + `SUPERNOVAE`, which the Makefile combines into
`STAR_FEEDBACK_ACTIVE` ([`Makefile:442`](../../Makefile)): a stellar-evolution lookup table
(mass-loss rates, wind velocities, SN ejecta, as a function of stellar mass and metallicity), read
by [`src/stars/star_tables.c`](../../src/stars/star_tables.c). The schema (from `star_tables.c`):
`Z_COUNT=13` (`Z_VALUES` from `1e-11` to `0.03`), `M_COUNT=32` (`M_VALUES` from `2` to `120` Msun),
and every `Z=.../M=.../` group carries `Age`, `Radius`, `Temperature`, `MassLossRate`, `HLossRate`,
`HeLossRate`, `MetalsLossRate`, `WindVelocity`, `SN_MassLoss`, `SN_HLoss`, `SN_HeLoss`,
`SN_MetalsLoss` (plus `Energy`/`Photons`, unused here since `RADIATION` is off). Units for these
fields are not documented in `src/stars/`; check with whoever built the table pipeline if feedback
strengths look physically off. The param files give it as `./star_tables/star_feedback_tables.hdf5`.

### `TreecoolFile` (COOLING)

Ships with the repo, no download needed: `TreecoolFile` points at `../../TREECOOL_ep` (an identical
copy is at [`data/TREECOOL_ep`](../../data/TREECOOL_ep)).

## 3. Build

You need a Grackle install (`GRACKLE_DIR` in `Makefile.systype`; the `Makefile` defaults to
`$(HOME)/software/grackle_solas/`) in addition to the usual MPI/GSL/HDF5 requirements in
[`documentation/source/running.md`](../../documentation/source/running.md).

```bash
make CONFIG=examples/agora_disc_bh_star_formation_3d/Config.sh \
     BUILD_DIR=examples/agora_disc_bh_star_formation_3d/build \
     EXEC=examples/agora_disc_bh_star_formation_3d/Arepo
```

## 4. Run

```bash
cd examples/agora_disc_bh_star_formation_3d
mkdir -p output_lowres                           # or output_medres
mpiexec -np <N> ./Arepo ./param_lowres.txt       # or ./param_medres.txt
```

BH quantities are in the `PartType5` block of the snapshots and in the `BLACKHOLES:` lines of the log.

## 5. Diagnostics

```bash
python agora_disc_diagnostics.py --snapshot output_lowres/snap_000.hdf5 output_lowres/snap_010.hdf5 \
    --output disc_profiles.png --grid-output "disc_grid_{name}.png"
```

Produces overlaid radial surface-density / scale-height / vertical velocity-dispersion profiles
across the given snapshots, plus one smoothed density-grid projection image per snapshot. Uses the
group's `analysistools` package if installed, otherwise a self-contained `h5py` reader; see the
script's docstring for all options (`--mask-radius`, `--ngrid`, `--save-data`, etc). For the
medium-res IC, remember the disc is PartType1 (section 1).

### Making the grid images smoother

The grid images deposit particles onto a 3D grid (CIC), smooth it with a Gaussian and plot the log
of the line-of-sight projection. If the result looks pixelated, speckled, or shows individual
particles or white holes, the grid is finer than the particle sampling. Tune:

| Option | Default | Effect |
|---|---|---|
| `--ngrid` | 256 | Cells per axis. Lower it (128-256) when particles are sparse. |
| `--smooth-sigma` | 3.0 | Gaussian width in cells. Should be at least the inter-particle spacing; raise to 4-6 for sparse regions. |
| `--interp` | `bilinear` | `imshow` interpolation (`nearest` for raw cells, `bicubic`/`gaussian` for smoother). |
| `--dex-range` | 5 | Colour scale spans this many dex below the maximum; empty cells are floored rather than drawn white. |
| `--projection` | `mean` | Line-of-sight reduction (`mean`, `sum`, `max`). `sum` gives a surface density. |
| `--cmap` | `viridis` | Colormap. |

### Gas images with vortrace

For gas (`--ptype 0`) you can skip the particle-to-grid step and ray-trace the simulation's own
Voronoi mesh with [vortrace](https://vortrace.readthedocs.io):

```bash
pip install vortrace
python agora_disc_diagnostics.py --snapshot output_lowres/snap_010.hdf5 --vortrace --ngrid 512
```

This gives a true column density, so it has no lattice or empty-cell artefacts and `--ngrid` can be
raised freely. It uses every gas cell, including background-grid cells, and needs `Density` and
`Masses` in `PartType0`. `--smooth-sigma` and `--projection` are ignored; `--interp`, `--dex-range`
and `--cmap` still apply. Other particle types and the radial profiles are unchanged. It is slower
than the grid deposit, and large background-grid cells can show up as streaks away from the disc.

A sparse low-resolution run without vortrace:

```bash
python agora_disc_diagnostics.py --snapshot output_lowres/snap_010.hdf5 \
    --ngrid 128 --smooth-sigma 3 --dex-range 4
```

## What has and has not been checked

Checked (on the previous versions of the ICs and parameters):

- The configuration builds (Open MPI, HDF5, GSL, Grackle) and the compile flags are as above.
- `get_ics.sh` produced valid files for both resolutions: header counts updated, BH ID above all
  others, BH mass 1e-3 code units, stellar-disc centre of mass within 0.02 kpc of the BH.

Not checked:

- **No simulation has been run.** An earlier attempt stopped at the missing star-feedback table, so
  IC reading, Bondi accretion and the feedback have not been exercised.
- **The new IC downloads and updated param files have not been run.** The new Dropbox URLs and the
  `get_ics.sh` changes (no unit-mass or `InternalEnergy` edits) are untested, and the BH position
  at `BoxSize/2` has not been rechecked against the new low-res IC.
- **Gas temperature.** Both ICs carry a gas `InternalEnergy` of ~101 (code units), used as supplied,
  so `InitGasTemp` in the param files has no effect. How Grackle treats this initial energy is
  untested.
- **Parameters:** the star-formation, feedback and metallicity values are the reference files'
  (see above), not validated against this fork's implementation. The BH parameters are
  copied from the cosmological examples or assumed.
- No `check.py`/reference data, unlike the `test.sh`-integrated examples: nothing here would catch a
  silent regression.

## Physics caveats to read before trusting results

- **Thermal feedback is continuous.** In `src/blackholes/bh_feedback.c` the energy is shared over
  all cells within the BH smoothing length, weighted by volume, every BH timestep, with no minimum
  heating temperature or stochastic coupling. In dense, cold gas the energy can be radiated away by
  Grackle cooling before it does any work, so the feedback may be weak at these resolutions.
- **Bondi rate is unboosted.** It uses the kernel-averaged density, sound speed and relative
  velocity within the BH radius, with no boost factor and no angular-momentum suppression. At this
  resolution the Bondi radius is far below the cell size, so the rate is a resolution-dependent
  estimate.
- **The BH is not pinned.** I found no repositioning in the BH code, so the BH moves as a normal
  collisionless particle. At low resolution it is lighter than a halo particle (1.25e7 Msun) and
  will be scattered; at medium resolution it is comparable (1.25e6 Msun). Expect some wandering
  from the centre, worst at low resolution.
- **BH softening** follows `SofteningTypeOfPartType5 = 1` (the halo softening), inherited from the
  parent.
- **No metal or mass return.** `Mload = 0` returns no mass to the gas.
