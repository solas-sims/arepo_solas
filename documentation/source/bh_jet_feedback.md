# Black hole jet feedback

`BH_JET_FEEDBACK` injects the kinetic power of a bipolar jet into the gas cells around a black hole. It is the
cell-based port of the virtual-particle jet used in Maragkakis et al. (2026); the jet is launched by the BH itself rather
than by spawned particles, so it works for any number of BHs, in cosmological runs, and with `BH_JET_PRESSURE_TRIGGER`.

Source: `src/blackholes/bh_jet_{common,density,feedback,trigger}.c`, the apply step in `bh_update.c`, and the timestep
limiter in `timestep.c`.

## How it fits together

Per BH step (`run.c`, `calculate_non_standard_physics_prior_mesh_construction()`), after `bh_density()`, `bh_accretion()` and `bh_swallow()`:

| Step | Routine | File | What it does |
|---|---|---|---|
| 1 | `bh_jet_pressure_trigger()` | `bh_jet_trigger.c` | Launch criterion: sets `BhP[].JetLaunchFlag` (only with `BH_JET_PRESSURE_TRIGGER`). |
| 2 | `bh_jet_prepare()` | `bh_jet_common.c` | Jet axis; adds this step's jet energy to the reservoir; decides how much to release. |
| 3a | `bh_jet_density()` + `bh_jet_feedback()` | `bh_jet_density.c`, `bh_jet_feedback.c` | Default cone mode: cone mass, then mass-weighted deposit into per-cell buffers. |
| 3b | `bh_jet_points()` + `bh_jet_points_inject()` | `bh_jet_points.c`, `bh_jet_points_inject.c` | `BH_JET_LAUNCH_POINT`: packets, flight delay, injection spheres, deposit into the same buffers. |
| 4 | `bh_jet_apply_cell()` | `bh_update.c` | At the end of the step, for each active cell with a buffer: applies the injection mode. |

Per-cell buffers are `SphP[].BhJetEnergy`, `BhJetMomentum[3]`, `BhJetMass` and the timestep mark `BhJetTiUntil`. Per-BH state, in
`BhP[]` (so it is written to restart files with the rest of the struct): `JetEnergyReservoir`, `JetEnergyRelease`, `JetConeMass`,
`JetAxis[3]`, and with `BH_JET_SPIN` `SpinVec[3]`, with `BH_JET_LAUNCH_POINT` `JetBurstCount` and the packet buffer `JetPacket[]`.
`bh_merger.c` adds the loser's reservoir (and packets' energy, and spin proxy) to the survivor. Buffers are applied the next time a
cell is active, so cells that will receive jet energy are marked for the timestep limiter first (`get_timestep_hydro()`).

## Configuration

| Option | Meaning |
| --- | --- |
| `BH_JET_FEEDBACK` | Enable the jet. Cannot be combined with `BH_THERMAL_FEEDBACK` (the accreted energy would be injected twice). |
| `MOMENTUM_JET`, `THERMAL_JET`, `ENERGY_JET` | Exactly one. See below. |
| `ENERGY_JET_EXACT` | With `ENERGY_JET`: the exact all-kinetic kick instead of the paper's eq. 12. |
| `BH_JET_FIXED_POWER` | Fixed power `BhJetPower` [erg/s] instead of `Epsilon_f Epsilon_r Mdot c^2`; creates the jet mass (arepobh behaviour). |
| `BH_JET_LAUNCH_POINT` | arepobh / paper geometry instead of the BH-centred cone: bursts fly out and are deposited in gas spheres at a launch radius (see below). |
| `BH_JET_SPIN` | Jet axis follows a spin proxy instead of the fixed axis. |
| `BH_JET_PRESSURE_TRIGGER` | Launch only if the jet ram pressure exceeds `BhJetPressureFactor` times the ambient pressure. Without it the jet always launches. |
| `BURST_MODE` | Release the stored energy only in multiples of `BhJetBurstMass * v_jet^2`. |
| `JET_TRACER` | Passive tracer set to 1 in every cell that receives jet material. |

Parameters: `BhJetVelocity` [km/s], `BhJetOpeningAngle` [deg, half-angle], `BhJetAxisX/Y/Z`, `BhJetTimestepFactor`,
`BhJetPressureFactor` (trigger), `BhJetLaunchRadius` [kpc], `BhJetDesNgb`, `BhJetDesDev` (launch-point mode), `BhJetPower` (fixed power), `BhJetBurstMass` [Msun] (burst mode), plus `Epsilon_f` and
`Epsilon_r` for the accretion-driven power.

## Model

**Power.** Per BH step, `E_new = Epsilon_f Epsilon_r Accretion c^2` (or `BhJetPower * dt` with `BH_JET_FIXED_POWER`).
Mass-loading `Mload` is not used by the jet.

**Reservoir.** `E_new` is added to `BhP[].JetEnergyReservoir`. If the launch condition holds the whole reservoir is
released (in `BURST_MODE`, the largest whole number of bursts). Otherwise it is carried over, so the time-averaged power
is not reduced by suppressed steps. The jet power on release can exceed the `L_jet` tested by the trigger, so the
trigger condition stays satisfied. Energy released into a cone that holds no gas returns to the reservoir.

**Cone.** Cells within the BH kernel radius `Hsml` and within `BhJetOpeningAngle` of the jet axis (either lobe) receive
the energy, in proportion to cell mass. The kick is radially outward from the BH, so the two lobes carry equal and
opposite momentum.

**Axis.** Fixed at `BhJetAxisX/Y/Z`, or with `BH_JET_SPIN` along `J`, where `J` starts as `M_BH * axis` and gains
`Accretion * L_gas / |L_gas|` each step (`L_gas` is the angular momentum of the kernel gas about the BH). This is a
proxy for the spin direction: it carries the history of the accreted gas's orientation, weighted by mass, not a spin
magnitude or a Bardeen-Petterson alignment.

**Injection** (applied the next time the cell is active):

* `MOMENTUM_JET`: add momentum `2E/v_jet`; thermal energy unchanged. Most of `E` is not deposited when the cell is
  more massive than the entrained jet mass `2E/v_jet^2`.
* `THERMAL_JET`: add the same momentum and add `E` to the total energy; what the momentum does not take up is heat.
* `ENERGY_JET`: add `E` to the total energy and kick along the jet direction by `sqrt(2 M K) - |p|`, with `K` the old kinetic
  energy plus `E` (Maragkakis et al. 2026, eq. 12, and arepobh). If `p` and the kick direction are not parallel the cell's kinetic
  energy ends up below `K`; the total energy has been raised by `E` regardless, so the shortfall becomes heat when the internal energy
  is recomputed. This is exact when the cell starts at rest or the kick is parallel to its momentum. `ENERGY_JET_EXACT` instead chooses
  the kick, `|p + s n|^2 / 2M = |p|^2 / 2M + E`, so that all of `E` is kinetic for any angle. In the example box the kinetic fraction of
  the energy gain was 0.630 with eq. 12 and 0.632 with the exact kick.

**Mass.** With `BH_JET_FIXED_POWER` the jet mass `2E/v_jet^2` is created in the cone cells, at rest and with no thermal
energy (as in arepobh). Otherwise the jet is made of entrained gas and cell masses are unchanged.

**Timestep.** Cells in the cone of a BH that holds jet energy get `dt <= BhJetTimestepFactor R_cell / (v_jet + c_s)`
for two BH steps after they were marked. Replaces the hard-coded radius-based limits of arepobh.

## Launch criterion (`BH_JET_PRESSURE_TRIGGER`)

A jet only launches if it can push its way out: its ram pressure at the launch base has to exceed the thermal pressure of the gas
it is launched into (`bh_jet_trigger.c`).

For each active BH, with `E_step` the jet energy of the step (accretion-driven: `Epsilon_f Epsilon_r Accretion c^2`, or
`BhJetPower dt` with `BH_JET_FIXED_POWER`) and `dt` its length, both in cgs:

* jet power `L_jet = E_step / dt`
* launch speed `v_jet = BhJetVelocity`; mass outflow rate `Mdot_j = 2 L_jet / v_jet^2` (a purely kinetic jet)
* launch area `A_jet = pi (R_inj sin theta)^2`, with `R_inj` the BH kernel radius `Hsml` (physical) and `theta = BhJetOpeningAngle`
* jet density `rho_jet = Mdot_j / (A_jet v_jet)`, so `P_ram = rho_jet v_jet^2 = 2 L_jet / (v_jet A_jet)`
* ambient pressure `P_amb`: the volume-weighted mean of the cell pressures in the BH kernel,
  `sum(P_cell V_cell) / sum(V_cell)`, accumulated in `bh_density()` (`BhP[].NgbsPressureVolume`, `NgbsVolume`)

`BhP[].JetLaunchFlag = 1` if `P_ram > BhJetPressureFactor * P_amb`. `BhJetPressureFactor = 1` is the plain criterion; larger values
make launching harder. `BhP[].JetRamPressure` and `AmbientPressure` hold the two pressures [cgs] for diagnostics (they are not in the
snapshots). The log reports, each BH step, how many active BHs can launch and the largest `P_ram / P_amb`:

```
BLACKHOLES: Jet pressure trigger: 1 of 1 active BHs can launch (max P_ram/P_amb = 5.4)
```

Units: everything is evaluated in cgs through the `cf_*` conversions (`cf_UnitTime_in_s`, `cf_UnitLength_in_cm`,
`cf_UnitPressure_in_cgs`), so it is also correct for comoving runs. The first version of the trigger took the BH timestep, which for
comoving runs is in ln(a), as a time and mixed code and physical units; that is fixed.

What the flag does depends on the mode:

* `BH_JET_FEEDBACK`: `bh_jet_prepare()` releases the jet energy only when the flag is 1. Energy produced while it is 0 is kept in
  `JetEnergyReservoir` and released at the next launch (in `BURST_MODE`, in whole bursts), so suppressed steps do not lower the mean
  power. On release the jet power is larger than the `L_jet` that was tested, so the condition is still satisfied.
* `BH_THERMAL_FEEDBACK` alone: the flag is only a diagnostic.
* Without `BH_JET_PRESSURE_TRIGGER`, `BH_JET_FEEDBACK` always launches.

The criterion always looks at the BH kernel, also in `BH_JET_LAUNCH_POINT` mode, where the gas that is actually hit lies at the
launch radius. Evaluating it at the injection sphere instead is not implemented.

## Launch-point geometry (`BH_JET_LAUNCH_POINT`)

This reproduces the injection model of Maragkakis et al. (2026), section 3.2, as implemented in arepobh. A released burst is a pair of virtual particles, carried as one *packet* (in `BURST_MODE` every whole burst is its own packet, with its own random direction, as in the paper's eq. 7):

1. **Direction.** Random within `BhJetOpeningAngle` of the jet axis, uniform in cos(theta). The draw is a function of the
   BH ID and the burst number only, so reruns and restarts give the same jets (arepobh used an unseeded `rand()` and a
   fixed 15 degrees; set `BhJetOpeningAngle 15` to match).
2. **Flight.** Each virtual particle flies in a straight line at `BhJetVelocity`, so it reaches the launch radius after
   `BhJetLaunchRadius / BhJetVelocity`. The packet stores the BH position at launch, and the flight is ballistic: the end point is
   the launch position plus R along +-d, whatever the BH does in the meantime. Up to
   `BH_JET_MAXPACKETS` (64) packets per BH are in flight; a full buffer folds new bursts into the newest packet.
3. **Delivery** (next BH step after arrival). Each lobe, at launch position +- R d, gets half the energy, shared in proportion to cell
   mass over the cells in a sphere whose radius holds `(BhJetDesNgb +- BhJetDesDev) * TargetGasMass` of gas (bisection, as in
   `bh_density`). The kick is along +-d. A sphere that cannot reach the target mass within R is
   accepted as is (arepobh stops the run). If either lobe's sphere is empty the packet is retried on the next step.
   The paper's fiducial runs (and arepobh `jets`) use the 10 nearest cells, a fixed number; the target-mass sphere is the paper's
   section 4.4 scheme (and arepobh `jets_improvements`). The paper finds that AREPO changes little between the two and PLUTO changes a lot.
4. **Timestep.** Cells at the target are marked when the packet is launched, so they already take short steps when it arrives.

The ambient-pressure trigger is unchanged: it still looks at the BH kernel, with the launch area at the kernel radius.
`JET_TRACER` cells (tracer above `BH_JET_TRACER_NODEREFINE`, 0.01) are never de-refined, as in arepobh `jets_improvements` (April 2026, after the paper; the paper only says that AREPO cells are kept within a factor of two of the target mass).

Tested in the same box with `BURST_MODE`, `BH_JET_FIXED_POWER`, refinement on: bursts launched at t = 1e-4, 2e-4, ... were
delivered 2e-4 later (the flight time for R = 2 kpc at 1e4 km/s); total energy gain equals the injected energy (90 = 9 bursts
x 10 to 1e-6); fast cells lie at r = 1.5-2.5 around R = 2; net momentum ~1e-16 of sum |p|. The de-refinement rule
and mesh splitting during a jet, are not tested. A restart with 12 packets in flight reproduces an uninterrupted run bitwise (4 tasks, `examples/bh_jet_box_3d`, repeated after the last changes to the packet record and the kick); a BH that migrates between tasks with packets in flight is not tested. With the BH given a velocity of 1e4 km/s across the box, the deposits stay at the launch position (not at the BH's current position): in that run the jet cells at the first delivery were centred 6 kpc behind the BH, at the launch x and at the launch radius in z. Whether the BH changed task during that run was not checked.

## Stratified atmosphere (`STATICBETA`)

The paper's stratified runs (section 4.5) use an isothermal beta model, `rho = rho_0 [1 + (r/r_c)^2]^(-3 beta/2)` (eq. 13), held in hydrostatic equilibrium by an
external acceleration `g = -3 beta k T / (mu m_p) r / (r^2 + r_c^2)` (eq. 14, with `(gamma-1) u = k T / (mu m_p)`). `STATICBETA` (in `grav_external.c`, with
`EXTERNALGRAVITY`) provides the acceleration and potential, centred on the box centre, with the constants in `Config.sh` (`BETA_MODEL_BETA`, `BETA_MODEL_RC` in code
length units, `BETA_MODEL_TEMP` in K, `BETA_MODEL_MU`). Non-comoving only. Unlike arepobh, it applies to all particle types, black holes included.
`examples/bh_jet_box_3d/create.py ... beta` makes the matching initial conditions (equal-mass cells sampled from the profile). In a 64 x 64 x 128 kpc box with
6.6e4 cells the gas kept its profile within 3 to 9 per cent in shells for 196 Myr, whereas without the acceleration the central density fell to a third in 49 Myr.
The 200 Myr relaxation of the paper and the remnant runs (jet switched off after 44 Myr) are not set up. See the example README.

## Differences from arepobh

* **Virtual particles are carried as records, not as particles.** In `BH_JET_LAUNCH_POINT` each burst is a pair of virtual particles
  (the +- lobes) that fly out and deposit their energy at the launch radius, as in the paper. arepobh implements them as
  type-5 particles spawned at the box centre, drifted at `JetVelocity` and destroyed when they pass 10 kpc. Here a burst is
  one *packet* (energy, direction, launch position, arrival time) in the launching BH's own `BhP[]` entry. The flight is straight
  and ballistic but is not followed step by step: the packet carries an arrival time of `BhJetLaunchRadius / BhJetVelocity`
  after launch, and the end point, launch position + R along +-d, is computed when it is delivered (arepobh launches from the fixed
  box centre). So no IDs, domain exchange or destruction are needed and any number of BHs work, but the virtual particles are not in
  snapshots, and delivery happens on the first BH step after the arrival time, not the step the particle crosses the radius. The default cone mode has no virtual particles: it deposits straight into a BH-centred cone with a radial kick.
* `ENERGY_JET` follows the paper's eq. 12 as arepobh does; the exact solve is the optional `ENERGY_JET_EXACT`.
* **Neighbour scheme.** Mass-targeted deposit sphere (paper section 4.4) rather than the 10 nearest cells (paper fiducial, arepobh `jets`). This was a choice, because the injection velocity then depends less on the local resolution.
* Dynamical timestep limiter; no hard-coded radii.
* Energy produced while the trigger is off is stored and released later.

## Checks

The test setup is a periodic 24^3 gas box with a BH at the centre (non-comoving, 2 MPI tasks):

* `ENERGY_JET`: the total energy gain equals the injected energy to the precision of `energy.txt` (~1e-5 relative),
  and after the injection stops it is constant; the kick gives exactly `sum(dKE) = sum(dE)`.
* Net momentum from the symmetric jet is ~1e-16 of `sum |p|`; all fast cells lie inside the cone.
* `THERMAL_JET`: total energy gain equals the injected energy. `MOMENTUM_JET`: gain equals the kinetic part only.
* A trigger that never fires injects nothing, and the reservoir equals `power * time`.
* `BURST_MODE`: releases are multiples of `BhJetBurstMass * v_jet^2`.
* `BH_JET_SPIN`: gas rotating about +x tilts the axis from +z towards +x when the accreted mass is made comparable to
  the BH mass (artificially, in a debug build).

Not tested: comoving integration (units follow the `cf_*` conversions used elsewhere), and restarts.
