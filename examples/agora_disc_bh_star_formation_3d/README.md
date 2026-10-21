# AGORA isolated disc + star formation + central black hole (3D)

The [`agora_disc_star_formation_3d`](../agora_disc_star_formation_3d/README.md) setup
(AGORA_SF star formation, Grackle chemistry, winds + supernovae) with a single central
supermassive black hole added:

- **BH in the ICs, no on-the-fly seeding.** A 1e6 Msun `PartType5` particle sits at the box
  centre (the AGORA discs are built about `BoxSize/2`), at rest in the construction frame.
  `BLACKHOLE_SEEDING`, `HALO_SEEDING` and FOF are *not* enabled.
- **Bondi accretion** (`BONDI_ACCRETION`), capped at the Eddington rate.
- **Thermal feedback only** (`BH_THERMAL_FEEDBACK`). Do not use `BH_FEEDBACK`: it also turns on the
  jet, which is not operational.
- Two resolutions, `param_lowres.txt` and `param_medres.txt`.

Like the parent example it is not wired into `test.sh`; run it by hand.

## What differs from the parent example

| | Parent | This example |
|---|---|---|
| `Config.sh` | -- | adds `BLACKHOLES`, `BONDI_ACCRETION`, `BH_THERMAL_FEEDBACK` |
| Parameters | -- | adds `BhDesNgb`, `BhDesDev`, `HMaxFactor`, `Epsilon_r`, `Epsilon_f`, `Mload` |
| IC | `..._without_u-with-grid.hdf5` | same file plus one `PartType5` (`add_central_bh.py`) |
| `param.txt` | one file | `param_lowres.txt`, `param_medres.txt` |

### Black hole parameters

| Parameter | Value | Meaning |
|---|---|---|
| `Epsilon_f` | 0.01 | fraction of the radiated luminosity coupled to the gas as thermal energy |
| `HMaxFactor` | 12.5 | maximum BH search radius = `HMaxFactor` x the BH's softening (type 5 -> softening type 1): 1.0 kpc low-res, 0.5 kpc med-res |
| `Epsilon_r` | 0.01 | radiative efficiency; sets the Eddington cap and the luminosity |
| `Mload` | 0 | mass-loading fraction. Energy deposited is `Epsilon_f (1 - Mload) Epsilon_r Mdot c^2`, mass returned `Mload Epsilon_r Mdot`, so 0 gives pure thermal feedback |
| `BhDesNgb`, `BhDesDev` | 32, 8 | smoothing-length targets, copied from the cosmological examples; untuned for a disc |

### Resolution-dependent values

| | low-res | med-res |
|---|---|---|
| gas / halo / disc / bulge particles | 1e5 / 1e5 / 1e5 / 1.25e4 | 1e6 / 1e6 / 1e6 / 1.25e5 |
| gas particle mass | 8.6e4 Msun | 8.6e3 Msun |
| `BoxSize` | 1200 kpc | 1523.4424 kpc (from the IC) |
| softenings (all types) | 0.08 kpc | 0.04 kpc (half) |
| `ReferenceGasPartMass` | 8.59322e-05 | 8.59322e-06 |
| `MeanVolume` | 52734.4 | 5273.44 |

The medium-resolution softenings are half the low-res ones, by choice. `ReferenceGasPartMass` and
`MeanVolume` are scaled by the particle-number ratio. None of these come from a reference medium-res run.
Everything else is as in the parent example.

## 1. Initial conditions

Requires `python3` with `h5py` and `numpy`.

### Low resolution

```bash
./get_ics.sh lowres
```

Downloads the parent example's low-res IC (already passed through `ADDBACKGROUNDGRID`) into
`./ICs/` and writes `./ICs/agora_lowres_arepo_ic_with_bh.hdf5`, with the gas `InternalEnergy` set to 0 so that
`InitGasTemp` (1e4 K) applies.

### Medium resolution

Uses the AGORA medium-resolution IC that has already been through `ADDBACKGROUNDGRID`
(`agora_galaxy_med_arepo-with-grid.hdf5`; not downloaded by the script). Put it in `./ICs/`, or pass its path:

```bash
./get_ics.sh medres [PATH_TO_IC]                 # writes ICs/agora_medres_arepo_ic_with_bh.hdf5
```

This file differs from the low-res one, and `get_ics.sh` adjusts the copy it writes:

| | medium-res file as supplied | in the written IC |
|---|---|---|
| header `UnitMass_in_g` | 1.989e43 (but its masses are in 1e9 Msun) | 1.989e42, matching `param_medres.txt` |
| gas `InternalEnergy` | disc gas ~2.5 (~245 K), grid cells 0 to 4.9e4 | 0, so `InitGasTemp` applies (as for low-res, below) |
| disc / halo particle types | PartType1 = stellar disc, PartType2 = halo | **unchanged** (low-res is the reverse) |
| `BoxSize` | 1523.44 (embedded parameters say 2000) | unchanged; `param_medres.txt` uses 1523.4424209594727 |
| BH position | -- | `BoxSize/2` = 761.72, which is the disc centre in this file |

The types are not swapped in the file. Softenings are equal for types 1 to 3, so the code is not affected, but
anything assuming PartType2 is the disc (for example `agora_disc_diagnostics.py`) must be pointed at PartType1.

`BH_MASS_MSUN=3e6 ./get_ics.sh ...` changes the BH mass.

## 2. Data tables, build, run

```bash
./get_tables.sh
make CONFIG=examples/agora_disc_bh_star_formation_3d/Config.sh \
     BUILD_DIR=examples/agora_disc_bh_star_formation_3d/build \
     EXEC=examples/agora_disc_bh_star_formation_3d/Arepo
cd examples/agora_disc_bh_star_formation_3d
mkdir -p output_lowres                           # or output_medres
mpiexec -np <N> ./Arepo ./param_lowres.txt       # or ./param_medres.txt
```

Diagnostics: `agora_disc_diagnostics.py`, as in the parent example. BH quantities are in the
`PartType5` block of the snapshots and in the `BLACKHOLES:` lines of the log.

## What has and has not been checked

Checked:

- The configuration builds (Open MPI, HDF5, GSL, Grackle) and the compile flags are as above.
- Both `param_lowres.txt` and `param_medres.txt` parse: every tag the code requires is present.
- `get_ics.sh` produces valid files for both resolutions: header counts updated, BH ID above all others, BH mass
  1e-3 code units, stellar-disc centre of mass within 0.02 kpc of the BH.

Not checked:

- **No simulation has been run.** The run stopped at the missing star-feedback table, so IC reading,
  Bondi accretion and the feedback have not been exercised.
- **Gas temperature.** The gas `InternalEnergy` is set to 0 in both written ICs, so `InitGasTemp = 10000`
  applies to every gas cell (neutral mean molecular weight, since the code only assumes full ionisation above
  1e4 K). The low-res file as downloaded has 2104.8735 in every massive gas cell (roughly 1e5 K by my
  arithmetic), which the parent example's README does not mention (it says the energy is 0). The first
  step has not been run, so how Grackle treats the 1e4 K start is untested.

## Physics caveats to read before trusting results

- **Thermal feedback is continuous.** In `src/blackholes/bh_feedback.c` the energy is shared over all
  cells within the BH smoothing length, weighted by volume, every BH timestep, with no minimum heating
  temperature or stochastic coupling. In dense, cold gas the energy can be radiated away by Grackle
  cooling before it does any work, so the feedback may be weak at these resolutions.
- **Bondi rate is unboosted.** It uses the kernel-averaged density, sound speed and relative velocity
  within the BH radius, with no boost factor and no angular-momentum suppression. At this resolution
  the Bondi radius is far below the cell size, so the rate is a resolution-dependent estimate.
- **The BH is not pinned.** I found no repositioning in the BH code, so the BH moves as a normal
  collisionless particle. At low resolution it is lighter than a halo particle (1.25e7 Msun) and will
  be scattered; at medium resolution it is comparable (1.25e6 Msun). Expect some wandering from the
  centre, worst at low resolution.
- **BH softening** follows `SofteningTypeOfPartType5 = 1` (the halo softening), inherited from the parent.
- **No metal or mass return.** `Mload = 0` returns no mass to the gas.
