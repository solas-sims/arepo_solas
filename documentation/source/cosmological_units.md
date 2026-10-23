# Cosmological Unit Conventions (comoving runs)

This page records how the scale factor `a` enters the solas modules, so that new physics code
converts between code and physical quantities consistently. With `ComovingIntegrationOn = 0`
every factor below is 1.

## Cosmology factors

`set_cosmo_factors_for_current_time()` (`src/time_integration/timestep.c`) refreshes these
whenever `All.Time` changes:

| Symbol | Value |
|---|---|
| `All.cf_atime` | `a` |
| `All.cf_a2inv`, `All.cf_a3inv` | `1/a^2`, `1/a^3` |
| `All.cf_hubble_a` | `H(a)` in code units |
| `All.cf_time_hubble_a` | `a H(a)` |

`All.Time` is the **scale factor** in comoving runs and the physical time otherwise. Never
subtract two `All.Time` values to get a duration; use `get_time_difference_in_Gyr(a0, a1)`.

## Storage conventions

| Quantity | Stored as | Physical value |
|---|---|---|
| Position, cell radius, `Hsml` | comoving code length | `L * cf_atime` |
| Gas density `SphP.Density` | comoving | `Density * cf_a3inv` |
| Velocity `P.Vel` | canonical, `a * v_pec` | `Vel / cf_atime` |
| Momentum | canonical, `a * p_phys` | `p / cf_atime` |
| Specific energy `Utherm`, sound speed | physical | as stored |
| Conserved `SphP.Energy` | `a^2 * E_phys` (thermal and kinetic) | `Energy * cf_a2inv` |
| Mass, halo mass | physical | as stored |
| Timestep (`dt` from time bins) | `d ln a` | physical `dt / cf_hubble_a` |

Rules that follow:

- Adding thermal energy to `SphP.Energy` needs `cf_atime^2 * du * Mass`. This applies to the SF
  temperature floor, cooling, and BH thermal feedback (`bh_update.c`).
- Free-fall times, Jeans masses and similar physical quantities use `Density * cf_a3inv`. A
  Jeans *length* compared against a cell radius is returned in comoving units,
  `sqrt(pi a / (G rho_comoving)) c_s`.
- Velocity differences combined with a sound speed (Bondi accretion) divide by `cf_atime` first.
- A physical timescale converts to the timestep variable with `dlna = H * dt`.

## `cf_Unit*` conversion constants (`set_cf_units()`)

Multiplying a code quantity by these gives physical CGS, including `a` and `h`:

- `cf_UnitLength_in_cm` carries `a/h`, `cf_UnitDensity_in_cgs` carries `h^2/a^3`.
- `cf_UnitMomentum_in_cgs` carries `1/a` and `cf_UnitEnergy_in_cgs` carries `1/a^2`, because the
  code stores canonical momentum and `a^2`-scaled energy. Convert physical feedback budgets with
  these (`code = phys / cf_Unit...`), not with `cf_UnitMass * cf_UnitVelocity`.
- `cf_UnitTime_in_s` / `cf_UnitTime_in_yr` convert a **`d ln a`** interval to physical time,
  since they carry `1/H`. They do not convert a difference of scale factors.

## Fixes applied (audit 2026-10-23)

The audit report with line references is `cosmo_expansion_factor_audit.md` in the repository
root. Fixed in this series:

- Free-fall time (`sfr_AGORA.c`, `sfr_JEANS.c`, `individual_star_formation.c`) now uses
  physical density.
- Jeans length and mass (`sfr_JEANS.c`) use the conventions above.
- Bondi / torque relative velocity (`bh_accretion.c`) is divided by `cf_atime`.
- BH thermal feedback energy (`bh_update.c`) gets the `cf_atime^2` factor.
- Stellar age (`star_update.c`) is the cosmic-time difference since birth, not `Delta a`.
- SN terminal-momentum cap (`star_feedback.c`) uses `cf_UnitMomentum_in_cgs`.
- SN lead time (`star_density.c`) is converted with `cf_UnitTime_in_yr`, the same units as
  `TimeToSN`.
- Bondi accretion (`bh_accretion.c`): kernel-summed gas density is converted to physical density
  (`* cf_a3inv`), and the rate, which is per code time, is divided by `cf_hubble_a` so that
  `rate * bh_timestep` (a `d ln a` interval) is a mass, as for the Eddington rate.
- `BH_CONSTANT_RADIUS` (`bh.h`, `BH_RADIUS_CODE`): `All.BhRadius` is now a *physical* radius and is
  converted to comoving code length (`/ cf_atime`) at every use.
- `star_radiation.c`: the RT ionisation and dissociation rates divide by `dt / cf_hubble_a` only.
  Grackle's `time_units = UnitTime/h` already absorbs `h`, so the extra `HubbleParam` was removed.
- `bh_feedback.c`: `h2` is now initialised and the unused `u = r * hinv` was removed (the file did
  not compile as written).

## Open items

- `bh_jet_feedback.c` has no cosmology handling, and the consuming injection is commented out
  (the code also references fields such as `KineticFeed` that are not in `SphP`). A note in the
  file lists what a re-enabled path needs: `cf_atime^2` on the energy and a canonical kick.
- These changes are untested in a comoving run with stars, SNe and black holes.
