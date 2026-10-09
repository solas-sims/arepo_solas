# Black hole jet feedback

`BH_JET_FEEDBACK` injects the kinetic power of a bipolar jet into the gas cells around a black hole. It is the
cell-based port of the virtual-particle jet used in Maragkakis et al. (2026); the jet is launched by the BH itself rather
than by spawned particles, so it works for any number of BHs, in cosmological runs, and with `BH_JET_PRESSURE_TRIGGER`.

Source: `src/blackholes/bh_jet_{common,density,feedback,trigger}.c`, the apply step in `bh_update.c`, and the timestep
limiter in `timestep.c`.

## Configuration

| Option | Meaning |
| --- | --- |
| `BH_JET_FEEDBACK` | Enable the jet. Cannot be combined with `BH_THERMAL_FEEDBACK` (the accreted energy would be injected twice). |
| `MOMENTUM_JET`, `THERMAL_JET`, `ENERGY_JET` | Exactly one. See below. |
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
* `ENERGY_JET`: add `E` to the total energy and kick along the cone direction by the momentum that puts all of `E` into
  kinetic energy, solving `|p + s n|^2 / 2M = |p|^2 / 2M + E` exactly for any angle between `p` and `n`. The arepobh
  form `|p'| - |p|` is only exact when the kick is parallel to the existing momentum.

**Mass.** With `BH_JET_FIXED_POWER` the jet mass `2E/v_jet^2` is created in the cone cells, at rest and with no thermal
energy (as in arepobh). Otherwise the jet is made of entrained gas and cell masses are unchanged.

**Timestep.** Cells in the cone of a BH that holds jet energy get `dt <= BhJetTimestepFactor R_cell / (v_jet + c_s)`
for two BH steps after they were marked. Replaces the hard-coded radius-based limits of arepobh.

## Launch-point geometry (`BH_JET_LAUNCH_POINT`)

This follows the `jets_improvements` branch of arepobh. A released burst is a *packet*:

1. **Direction.** Random within `BhJetOpeningAngle` of the jet axis, uniform in cos(theta). The draw is a function of the
   BH ID and the burst number only, so reruns and restarts give the same jets (arepobh used an unseeded `rand()` and a
   fixed 15 degrees; set `BhJetOpeningAngle 15` to match).
2. **Flight.** The packet takes `BhJetLaunchRadius / BhJetVelocity` to reach the launch radius. Up to
   `BH_JET_MAXPACKETS` (64) packets per BH are in flight; a full buffer folds new bursts into the newest packet.
3. **Delivery** (next BH step after arrival). Each lobe, at BH position +- R d, gets half the energy, spread over a
   sphere whose radius holds `(BhJetDesNgb +- BhJetDesDev) * TargetGasMass` of gas (bisection, as in `bh_density`),
   in proportion to cell mass, with the kick along +-d. A sphere that cannot reach the target mass within R is
   accepted as is (arepobh stops the run). If either lobe's sphere is empty the packet is retried on the next step.
4. **Timestep.** Cells at the target are marked when the packet is launched, so they already take short steps when it arrives.

The ambient-pressure trigger is unchanged: it still looks at the BH kernel, with the launch area at the kernel radius.
`JET_TRACER` cells (tracer above `BH_JET_TRACER_NODEREFINE`, 0.01) are never de-refined, as in arepobh.

Tested in the same box with `BURST_MODE`, `BH_JET_FIXED_POWER`, refinement on: bursts launched at t = 1e-4, 2e-4, ... were
delivered 2e-4 later (the flight time for R = 2 kpc at 1e4 km/s); total energy gain equals the injected energy (90 = 9 bursts
x 10 to 1e-6); fast cells lie at r = 1.5-2.5 around R = 2; net momentum ~1e-16 of sum |p|. The de-refinement rule
and mesh splitting during a jet, and restarts with packets in flight, are not tested.

## Differences from arepobh

* Jets are launched by BHs, not from virtual particles at the box centre. The default cone mode kicks radially in a BH-centred cone; `BH_JET_LAUNCH_POINT` reproduces the arepobh deposition geometry.
* `ENERGY_JET` uses the exact momentum solve above.
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
