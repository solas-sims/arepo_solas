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
* **Gas:** uniform and at rest, rho_0 = 2.5e-26 g/cm^3 and T = 10^7 K, the uniform medium of Maragkakis et al. (2026), section 3.3,
  on a 32 x 32 x 64 grid (dx = 2 kpc; the paper uses a regular grid of 10^8 elements). There is no gravity, so a uniform medium
  is the only equilibrium state. The paper's stratified atmosphere is in `examples/bh_jet_stratified_3d`.
* **Units:** 1 kpc, 1 km/s, 1e6 Msun (so the time unit is 0.978 Gyr), as in the arepobh `param.txt`.
* **Jet:** 1e45 erg/s at 1e4 km/s, in bursts of `BhJetBurstMass * v_jet^2` = 2e57 erg (every ~63 kyr), within 15 degrees of
  the z axis. A burst takes `BhJetLaunchRadius / v_jet` (0.5 Myr) to reach 5 kpc and is then spread over the gas in a
  sphere around each lobe point that holds 8 +- 2 target gas masses (the neighbour-mass scheme, see below). `ENERGY_JET` follows the paper's eq. 12: the kick is
  sqrt(2 M K) - |p| along the direction of travel, and any shortfall in kinetic energy becomes heat.
* **BH:** 1e6 Msun. `BONDI_ACCRETION` is on only because `BH_FEEDBACK_ACTIVE` needs an accretion model; with
  `BH_JET_FIXED_POWER` the jet does not use the accretion rate.
* **Run length:** `TimeMax 0.003` is 2.9 Myr. The BH's timestep is limited by its neighbours', which are long until the
  first burst lands (7.5e-4 here), so the jet energy is produced in three chunks (one per BH step) and 34 bursts
  are released, not the 46 that 1e45 erg/s over 2.9 Myr would give. 29 are delivered by the end and 5 are
  still in flight. Each burst is its own pair of virtual particles with its own random direction. Cells in an active jet region take timesteps down to ~2e-5.

## 1. Initial conditions

```bash
mkdir -p run && cp Config.sh param.txt run/
python3 create.py run 32
```

The second argument is the number of cells per side in x and y (default 32; z gets twice as many). It prints the
`ReferenceGasPartMass` that matches the grid (2.96 for 32); change it in `param.txt` if you change the resolution.
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
   code reports as taken up by the gas, to 1e-3. In the reference run the difference is 1e-13.
2. **Momentum.** The net momentum is below 1% of the sum of |p|. The two lobes are launched symmetrically, but the
   kick depends on each cell's own momentum, so it is not exactly zero (2e-3 in the reference run).
3. **Tracer.** Cells carry jet tracer, and they lie around the launch radius and beyond.

It also writes `jet_slice.png`: |v| and the jet tracer in a slab through the BH, showing the two lobes along z.

Reference run (4 tasks, `-O2`): 29 bursts' worth of energy (2.9e9 code units) delivered by t = 0.003, 265 cells with
tracer at 2-13 kpc from the BH. The jet energy in flight at the end (5 bursts) is about the power times the 0.5 Myr flight time, as expected.

The energy bookkeeping lines to look at in `log.txt`:

```
BLACKHOLES: Jet energy [code]: produced ..., released ... by 1 of 1 active BHs, held in reservoirs ...
BLACKHOLES: Jet launch points: N packets delivered, ... postponed (empty sphere), ... spheres capped at the launch radius, ... in flight
BLACKHOLES: Energy given by BHs = ..., Energy taken up by gas particles = ...
```

`given` and `taken` should be equal; `postponed` and `capped` should be 0 here.

## Relation to the paper runs

The paper is Maragkakis et al. (2026, MNRAS 547, stag368). Its fiducial jet runs correspond to the arepobh `jets` branch (July 2025). The
`jets_improvements` branch dates from April 2026, after the paper was accepted. The settings below are from the `jets` files and agree with the paper's
section 3.3; the example differs in the last column.

| | arepobh `param.txt` / `Config.sh` | This example |
|---|---|---|
| Box | 300 kpc, `LONG_Z=2` (300 x 300 x 600), periodic hydro | 64 x 64 x 128 kpc |
| Duration | `TimeMax 0.1` (98 Myr) | 0.003 (2.9 Myr) |
| Jet power | `JetPower 5e45` erg/s | 1e45 |
| Jet speed | `JetVelocity 4e4` km/s | 1e4 |
| Burst mass | `JetParticleMass 0.04` (4e4 Msun) | `BhJetBurstMass 1e6` Msun |
| Launch radius | 10 kpc (hard-coded in `bh_virtual_particles.c`) | `BhJetLaunchRadius 5` |
| Cone | 15 degrees (hard-coded) | `BhJetOpeningAngle 15` |
| Deposit region | the 10 nearest gas elements (`BhDesNgbMass 10`, which counts cells, +- 2) | a sphere holding `BhJetDesNgb` 8 +- 2 target gas masses (the paper's section 4.4 neighbour-mass scheme) |
| Resolution | `ReferenceGasPartMass 0.2` (2e5 Msun, N = 10^8) | 2.96 (3e6 Msun) |
| Flags | `ENERGY_JET`, `BURST_MODE`, `LONG_Z=2`, mass-based refinement | same, plus `BH_JET_FIXED_POWER`, `BH_JET_LAUNCH_POINT`, `JET_TRACER` |

To get closer to the paper, set `BoxSize 300`, `TimeMax 0.1`, `BhJetPower 5e45`, `BhJetVelocity 4e4`, `BhJetBurstMass 4e4`,
`BhJetLaunchRadius 10`, `BhJetDesNgb 10`, `ReferenceGasPartMass 0.2` and make a larger grid with `create.py` (its gas is already the paper's
uniform medium). The paper's fiducial run has 10^8 elements, so this needs a cluster. Checks against the paper's uniform-medium results (its Figs 5-7, AREPO curves at 98 Myr): lobe length about 250 kpc, lobe aspect ratio about 4.8,
and about 60 per cent of the injected energy in the ambient medium, mostly as heat. Lobe cells are those with T > 10^8 K, v > 5000 km/s, or tracer > 0.01.

## Variations

In `Config.sh`:

* Remove `BH_JET_LAUNCH_POINT` for the BH-centred cone: the energy then goes straight into the cone cells inside the BH
  kernel (see `documentation/source/bh_jet_feedback.md`).
* Replace `ENERGY_JET` by `THERMAL_JET` or `MOMENTUM_JET` for the other two injection modes (exactly one must be set).
* Remove `BH_JET_FIXED_POWER` to drive the jet from `Epsilon_f * Epsilon_r * Mdot * c^2` (drop `BhJetPower` from
  `param.txt`; the Bondi rate for a 1e6 Msun BH in this gas is tiny, so give it a much heavier BH or a denser medium).
* Add `BH_JET_STOP_TIME` (and `BhJetStopTime`, code time) to switch the jet off at a given time; see `examples/bh_jet_stratified_3d`.
* Add `BH_JET_PRESSURE_TRIGGER` (and `BhJetPressureFactor` to `param.txt`) to launch only when the jet ram pressure
  exceeds the ambient pressure; `BH_JET_SPIN` to let the jet axis follow the accreted angular momentum.

## Running on many tasks: `OVERRIDE_PEANOGRID_WARNING`

`Config.sh` sets `OVERRIDE_PEANOGRID_WARNING`. Without it, a run on many tasks (32 in the test) stops with

```
DOMAIN: Note: we would like to refine top-tree, but PEANOGRID is not fine enough
TERMINATE: ... domain_determineTopTree() ... Consider setting BITS_PER_DIMENSION up to a value of 42 ...
```

`BITS_PER_DIMENSION` already defaults to 42, so the hint does not apply. The cause is the domain workload: the top-tree must split every
node whose share of the work exceeds `1 / (TopNodeFactor * MultipleDomains * NTask)` (1/256 for 32 tasks here, 1/32 for 4). Cells in the jet
region are in the shortest time bins and carry far more work than the others, so at a sync-point where only a few cells are active a single
cell can exceed that share, and it cannot be split. The override lets the run continue; the cost is one less well-balanced domain
decomposition. In the 32-task test the note printed once, the run finished, and `check.py` gave the same results as on 4 tasks (energy
closure 3e-13, 206 against 208 tracer cells, since the decomposition differs; measured before the later changes to the kick formula and the burst splitting). This is not a coordinate problem (NaN or coincident particles), which would be
a different failure to look for.

## Restart test

Packets in flight are stored in the BH's `BhP` entry, so they go into the restart files. To check this, stop the run while a
packet is in flight and restart it:

```bash
cd run
mpirun -np 4 ./Arepo param.txt > log1.txt &     # start the run
touch output/stop                                # once the log shows Time > 0.00165 and some packets "in flight"
# (wait for it to exit: it writes restart files and stops)
mpirun -np 4 ./Arepo param.txt 1 > log2.txt      # restart
```

In the reference test (stopped at t = 0.00171 with 12 packets in flight, restarted with the same 4 tasks) the restarted run's
launch and delivery events were identical to an uninterrupted run, and all five snapshots, including the one after the restart, were
bitwise identical (positions, velocities, masses, internal energies and tracer). The packets in flight at the stop were delivered after the
restart.

## Not tested here

A restart on a different number of tasks (Arepo needs the same number as the run that wrote the restart files), a BH that moves
between tasks while packets are in flight (the BH sits still in this example), and comoving integration.
