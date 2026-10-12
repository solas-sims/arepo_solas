# BH jet in a hot-gas box (3D)

A smoke test and starting point for `BH_JET_FEEDBACK` with the `BH_JET_LAUNCH_POINT` geometry, which follows the
arepobh runs of Maragkakis et al. (2026). A black hole sits at the centre of a uniform, hot, periodic box and
launches bursts that are deposited in gas spheres 5 kpc up and down the z axis. The set-up is scaled down so it
runs in a couple of minutes; the paper's values are listed below.

It does **not** plug into `test.sh`; run it by hand as described below.

## What is in this directory

| File | Purpose |
|---|---|
| `Config.sh` | Compile-time flags |
| `param.txt` | Runtime parameters |
| `create.py` | Writes `IC.hdf5` (needs `numpy` and `h5py`) |
| `check.py` | Post-run checks (energy accounting, momentum, tracer); optional plot |

## The set-up

* **Box:** `BoxSize 64` with `LONG_Z=2`, so 64 x 64 x 128 kpc. `BoxSize` is the x and y length; z is `LONG_Z` times
  that. The BH is at the centre, (32, 32, 64) kpc, and the jet axis is z.
* **Boundaries:** the **hydro box is periodic**. The `REFLECTIVE_*` flags are not set, exactly as in the arepobh
  `Config.sh`. `GRAVITY_NOT_PERIODIC` is set and `param.txt` has `PeriodicBoundariesOn 0`, but there is no
  `SELFGRAVITY` and `GravityConstantInternal 0`, so this only satisfies the consistency check in
  `check_parameters()`. Do not read `PeriodicBoundariesOn 0` as open boundaries for the gas.
* **Gas:** uniform and at rest, n_H = 0.01 cm^-3, T = 10^7 K, on a 32 x 32 x 64 grid (dx = 2 kpc). There is no
  gravity, so a uniform medium is the only equilibrium state. This is not the paper's initial condition (see below).
* **Units:** 1 kpc, 1 km/s, 1e6 Msun (so the time unit is 0.978 Gyr), as in the arepobh `param.txt`.
* **Jet:** 1e45 erg/s at 1e4 km/s, in bursts of `BhJetBurstMass * v_jet^2` = 2e57 erg (every ~63 kyr), within 15 degrees of
  the z axis. A burst takes `BhJetLaunchRadius / v_jet` (0.5 Myr) to reach 5 kpc and is then spread over the gas in a
  sphere around each lobe point that holds 8 +- 2 target gas masses. `ENERGY_JET` makes all of it kinetic energy.
* **BH:** 1e6 Msun. `BONDI_ACCRETION` is on only because `BH_FEEDBACK_ACTIVE` needs an accretion model; with
  `BH_JET_FIXED_POWER` the jet does not use the accretion rate.
* **Run length:** `TimeMax 0.003` is 2.9 Myr. The BH's timestep is limited by its neighbours', which are long until the
  first burst lands (7.5e-4 here), so the jet energy is produced in three chunks (one per BH step) and about 35 bursts
  are released, not the 46 that 1e45 erg/s over 2.9 Myr would give. About 27 are delivered by the end and 7-8 are
  still in flight. Cells in an active jet region take timesteps down to ~2e-5.

## 1. Initial conditions

```bash
mkdir -p run && cp Config.sh param.txt run/
python3 create.py run 32
```

The second argument is the number of cells per side in x and y (default 32; z gets twice as many). It prints the
`ReferenceGasPartMass` that matches the grid (2.6 for 32); change it in `param.txt` if you change the resolution.
`BhJetDesNgb` is in units of that mass.

## 2. Build

From the repository root:

```bash
make -j8 SYSTYPE=MACOSX CONFIG=examples/bh_jet_box_3d/run/Config.sh BUILD_DIR=examples/bh_jet_box_3d/run/build EXEC=examples/bh_jet_box_3d/run/Arepo
```

(`SYSTYPE` as for your machine.) The default macOS/Linux flags are `-O0 -g`, which makes this run about an order of
magnitude slower. For a quick run add `OPTIMIZE="-std=c11 -O2 -g -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function"`
to the `make` line.

## 3. Run

```bash
cd run
mpirun -np 4 ./Arepo param.txt > log.txt
```

With `-O2` and 4 tasks the run takes about 1.5 minutes and writes five snapshots. Keep the stdout in `log.txt`:
`check.py` reads the jet energy accounting from it.

## 4. Check

```bash
python3 check.py run
```

It compares the last snapshot with the IC and exits non-zero if a check fails:

1. **Energy closes.** The total energy gain (from the snapshots, at full precision) equals the jet energy the
   code reports as taken up by the gas, to 1e-3. In the reference run the difference is 3e-13.
2. **Momentum.** The net momentum is below 1% of the sum of |p|. The two lobes are launched symmetrically, but the
   kick depends on each cell's own momentum, so it is not exactly zero (7e-4 in the reference run).
3. **Tracer.** Cells carry jet tracer, and they lie around the launch radius and beyond.

It also writes `jet_slice.png`: |v| and the jet tracer in a slab through the BH, showing the two lobes along z.

Reference run (2026-10-12, 4 tasks, `-O2`): 27 bursts' worth of energy (2.7e9 code units) delivered by t = 0.003, 208 cells with
tracer at 4-14 kpc from the BH. The jet energy in flight at the end (~0.8e9) is the power times the 0.5 Myr flight time, as expected.

The energy bookkeeping lines to look at in `log.txt`:

```
BLACKHOLES: Jet energy [code]: produced ..., released ... by 1 of 1 active BHs, held in reservoirs ...
BLACKHOLES: Jet launch points: N packets delivered, ... postponed (empty sphere), ... spheres capped at the launch radius, ... in flight
BLACKHOLES: Energy given by BHs = ..., Energy taken up by gas particles = ...
```

`given` and `taken` should be equal; `postponed` and `capped` should be 0 here.

## Relation to the paper runs

The arepobh `jets_improvements` branch (the version used for the paper) has these settings. The example differs in the
last column.

| | arepobh `param.txt` / `Config.sh` | This example |
|---|---|---|
| Box | 300 kpc, `LONG_Z=2` (300 x 300 x 600), periodic hydro | 64 x 64 x 128 kpc |
| Duration | `TimeMax 0.1` (98 Myr) | 0.003 (2.9 Myr) |
| Jet power | `JetPower 5e45` erg/s | 1e45 |
| Jet speed | `JetVelocity 4e4` km/s | 1e4 |
| Burst mass | `JetParticleMass 0.04` (4e4 Msun) | `BhJetBurstMass 1e6` Msun |
| Launch radius | 10 kpc (hard-coded in `bh_virtual_particles.c`) | `BhJetLaunchRadius 5` |
| Cone | 15 degrees (hard-coded) | `BhJetOpeningAngle 15` |
| Sphere gas mass | `BhDesNgbMass 10` target masses +- 2 | `BhJetDesNgb 8` +- 2 |
| Resolution | `ReferenceGasPartMass 0.2` (2e5 Msun) | 2.6 (2.6e6 Msun) |
| Flags | `ENERGY_JET`, `BURST_MODE`, `LONG_Z=2`, mass-based refinement | same, plus `BH_JET_FIXED_POWER`, `BH_JET_LAUNCH_POINT`, `JET_TRACER` |

To get closer to the paper, set `BoxSize 300`, `TimeMax 0.1`, `BhJetPower 5e45`, `BhJetVelocity 4e4`, `BhJetBurstMass 4e4`,
`BhJetLaunchRadius 10`, `BhJetDesNgb 10`, `ReferenceGasPartMass 0.2` and make a larger grid with `create.py`. That is a much
larger run (~54 million cells at the paper's resolution if uniform; the paper used refinement) and needs a cluster. The
paper's initial conditions are not in the arepobh repository; whether its gas was uniform or a profile is not known from the
files, so this example does not claim to reproduce it.

## Variations

In `Config.sh`:

* Remove `BH_JET_LAUNCH_POINT` for the BH-centred cone: the energy then goes straight into the cone cells inside the BH
  kernel (see `documentation/source/bh_jet_feedback.md`).
* Replace `ENERGY_JET` by `THERMAL_JET` or `MOMENTUM_JET` for the other two injection modes (exactly one must be set).
* Remove `BH_JET_FIXED_POWER` to drive the jet from `Epsilon_f * Epsilon_r * Mdot * c^2` (drop `BhJetPower` from
  `param.txt`; the Bondi rate for a 1e6 Msun BH in this gas is tiny, so give it a much heavier BH or a denser medium).
* Add `BH_JET_PRESSURE_TRIGGER` (and `BhJetPressureFactor` to `param.txt`) to launch only when the jet ram pressure
  exceeds the ambient pressure; `BH_JET_SPIN` to let the jet axis follow the accreted angular momentum.

## Not tested here

Restarts with jet packets in flight, and comoving integration.
