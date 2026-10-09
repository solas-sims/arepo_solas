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
| `BH_JET_SPIN` | Jet axis follows a spin proxy instead of the fixed axis. |
| `BH_JET_PRESSURE_TRIGGER` | Launch only if the jet ram pressure exceeds `BhJetPressureFactor` times the ambient pressure. Without it the jet always launches. |
| `BURST_MODE` | Release the stored energy only in multiples of `BhJetBurstMass * v_jet^2`. |
| `JET_TRACER` | Passive tracer set to 1 in every cell that receives jet material. |

Parameters: `BhJetVelocity` [km/s], `BhJetOpeningAngle` [deg, half-angle], `BhJetAxisX/Y/Z`, `BhJetTimestepFactor`,
`BhJetPressureFactor` (trigger), `BhJetPower` (fixed power), `BhJetBurstMass` [Msun] (burst mode), plus `Epsilon_f` and
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

## Differences from arepobh

* Jets are launched by BHs, not from virtual particles at the box centre; the cone direction is the BH axis, not random.
* Radial kick inside the cone (arepobh kicked along the axis).
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
