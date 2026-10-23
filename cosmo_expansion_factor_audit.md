# Cosmological Expansion-Factor Audit — arepo_solas-fof

Audit date: 2026-09-08
Branch: `merge-fof-into-star-feedback`
Scope: where the scale factor `a` (`All.Time`, `All.cf_atime`, `All.cf_a2inv`, `All.cf_a3inv`, `All.cf_hubble_a`, etc.) is required for correct comoving <-> physical unit conversion, and where it looks missing. Audit only — no files were modified.

---

## 1. Baseline convention (confirmed correct, established elsewhere in the code)

Setup: `set_cosmo_factors_for_current_time()` in `src/time_integration/timestep.c:77-108`, called whenever `All.Time` changes (from `begrun.c:205`, `init.c:106`, `read_ic.c:990`, `io.c:609`).

```
All.cf_atime         = a                       (timestep.c:81)
All.cf_a2inv         = 1/a^2                    (timestep.c:82)
All.cf_a3inv         = 1/a^3                    (timestep.c:83)
All.cf_afac1/2/3     = pow(a, ...GAMMA...)      (timestep.c:84-86)   -- hydro EOS/pressure factors
All.cf_hubble_a      = H(a)                     (timestep.c:87)
All.cf_time_hubble_a = a*H(a)                   (timestep.c:88)
```

Fork-specific addition, `set_cf_units()` (`timestep.c:110-140`), bakes `a` (and `h`) scaling directly into `cf_UnitDensity_in_cgs` (∝1/a³), `cf_UnitMomentum_in_cgs` (∝1/a), `cf_UnitEnergy_in_cgs` (∝1/a²).

**Per-quantity convention** (cross-checked across many call sites):

| Quantity | Storage convention | Physical value | Verified correct at |
|---|---|---|---|
| Density (`SphP.Density`) | comoving | `Density * All.cf_a3inv` | `sfr_eEOS.c:118,121,125,138,226,229,232,237`; `cooling.c:929` |
| Velocity (`P.Vel`) | canonical (`Vel = a * v_peculiar`) | `Vel / All.cf_atime` | `timestep.c:427,464-465`; `subfind_so.c:729,777,819` (Hubble-flow-corrected: `Vel[i]/cf_atime - Vel_centre + cf_Hrate*Pos_pbc`) |
| Position/length | comoving code length | `length * All.cf_atime` | `subfind_so.c:769-771,811-813`; `update_primitive_variables.c:144` |
| Internal energy (Utherm) | physical (a-independent) already | n/a | consistent throughout |
| Conserved Energy increment | needs `cf_atime^2` factor when incrementing from a Utherm change | — | `sfr_eEOS.c:114,179`; `sfr_AGORA.c:126`; `sfr_JEANS.c:164`; `cooling.c:941`; `sfr_starbystar.c:109` |
| Potential/acceleration | standard comoving-potential / Hubble-flow terms | — | `accel.c:214-343` (unmodified upstream Arepo, not touched by recent commits) |
| Overdensity threshold (`All.OverDensThresh`) | compared against raw **comoving** Density (deliberate — mean baryon density is const. in comoving units), guarded by `if(All.ComovingIntegrationOn)` | — | `starformation.c:83-86,142-144`; `sfr_eEOS.c:142-144,236-238`; `sfr_AGORA.c:52-54`; `sfr_JEANS.c:91-93` |

---

## 2. Checked and confirmed correct (no action needed)

- `src/blackholes/bh_seed.c:184-185` and `src/fof/fof.c:832-875,1008-1032` — halo/BH seed velocity is a mass-weighted average of `P[].Vel`, all in the same canonical convention; no Hubble-flow term needed since it isn't a spatially-offset relative velocity.
- `sf_threshold_halo_mass_factor()` (`starformation.c:122-128`) and `SphP.HostHaloMass` / `All.MinHaloMassForNormalSF` — halo mass is a-independent; comparison is unit-consistent.
- `sfr_AGORA.c:52-53,62-64`, `sfr_JEANS.c:91-92` — comoving-overdensity floor and physical-density conversion via `All.cf_UnitDensity_in_cgs` both correctly scaled.
- `src/blackholes/bh_accretion.c:344-346` (position difference in Bondi kernel loop) — no `cf_atime` applied, but acceptable since it's only used inside the SPH kernel radius for weighting, not compared against an absolute physical length.

---

## 3. Suspected / confirmed missing expansion-factor conversions

### (a) Freefall time computed from raw comoving Density instead of physical density
Same bug, three independent locations:

- `src/star_formation/sfr_AGORA.c:160`
- `src/star_formation/sfr_JEANS.c:198`
- `src/star_formation/individual_star_formation/individual_star_formation.c:94`

All compute:
```c
t_freefall = sqrt(3*M_PI/32/All.G/SphP[i].Density)
```
without multiplying `SphP[i].Density` by `All.cf_a3inv`. Contrast with the same file's own `sf_criteria()` (`sfr_AGORA.c:62`), which *does* convert density to physical via `cf_UnitDensity_in_cgs`. Also contrast with the reference-correct pattern in `sfr_eEOS.c:121/229`, where `calc_egyeff(i, dens * All.cf_a3inv, ...)` feeds physical density into an equivalent `tsfr = sqrt(PhysDensThresh/rho)*MaxSfrTimescale`.

**Effect**: `t_freefall` (and hence SFR normalization) is off by `a^{3/2}` — wrong at high z (a<1).

### (b) Jeans length and Jeans mass — same missing `cf_a3inv`, but feeds the SF criterion directly, not just a rate
- `src/star_formation/sfr_JEANS.c:62` (`get_jeans_length`)
- `src/star_formation/sfr_JEANS.c:79` (`get_jeans_mass`)

Both use raw `SphP[i].Density` inside `sqrt(.../Density)`. These are compared at:
- `sfr_JEANS.c:107` — `jeans_mass < threshold*Mass` (against `P[i].Mass`, genuinely physical/a-independent)
- `sfr_JEANS.c:116` — `jeans_length < 2*cell_radius` (against a genuinely comoving length, per `update_primitive_variables.c:144`'s explicit `cf_atime` multiplication)

So both comparisons currently mix scale-factor conventions. Derived correction:
- correct comoving Jeans length = `sqrt(a) *` (current value)
- correct physical Jeans mass = (current value) `* a^{1.5}`

As written, both are over-estimated at high z (a<1), making the `JEANS_SF` gate harder to satisfy than intended early in a run — this is a **criterion-level** bug, not merely a normalization bug.

### (c) BH Bondi accretion relative velocity not converted canonical -> peculiar
- `src/blackholes/bh_accretion.c:351-353`:
```c
dvx = P[i].Vel[0] - vel[0];   // and dvy, dvz
```
computed from raw canonical `Vel`, no division by `All.cf_atime`. Feeds `gas_velocity[]` (density-weighted accumulator, `bh_accretion.c:379-381`) -> `BhP[i].GasVelocity` -> `gas_velocity_norm` at `bh_accretion.c:489-490`, combined with `gas_sound_speed` (correctly scale-independent, since pressure/density ratio cancels the comoving normalization) in the Bondi denominator at `bh_accretion.c:492` (`c_s^2 + v^2`).

Mixing a canonical velocity term with a physical sound-speed term in the same sum is inconsistent. Established fix pattern: divide by `All.cf_atime` before use (per `subfind_so.c:777,819`, `timestep.c:427,464-465`).

**Effect**: at a<1 (high z), the velocity term in the Bondi denominator is under-estimated relative to its true physical value, systematically biasing `BondiRate` upward at high z whenever the velocity term isn't negligible relative to sound speed.

Same missing `/All.cf_atime` propagates into `TORQUE_ACCRETION`'s `v_cross`/`v_phi` (`bh_accretion.c:368-371,393`), though there the effect is only on a sign/threshold test (`is_disk = v_phi>0`), not a magnitude — lower severity, worth a quick look but not urgent.

### (d) BH thermal-feedback energy injected without the `cf_atime^2` conserved-energy factor
- `src/blackholes/bh_feedback.c:177,226`: computes a **physical** energy
```c
energy_feed = Epsilon_f*(1-Mload)*Epsilon_r*accretion*(CLIGHT^2/cf_UnitVelocity_in_cm_per_s^2)
```
(`cf_UnitVelocity_in_cm_per_s` itself has no `a`-dependence, per `timestep.c:125`). Accumulated into `SphP[i].BhThermalFeed`.

- `src/blackholes/bh_update.c:287`:
```c
SphP[i].Energy += SphP[i].BhThermalFeed;
```
adds this physical energy **directly** with no `All.cf_atime^2` factor, unlike every other thermal-energy-injection site in the codebase (`sfr_eEOS.c:114,179`; `sfr_AGORA.c:126`; `sfr_JEANS.c:164`; `cooling.c:941`; `sfr_starbystar.c:109`), all of which multiply by `cf_atime^2` before adding to `Energy`.

**Effect**: if genuine (not an intentional convention choice), BH thermal feedback is systematically too weak at high z (a<1) relative to the accretion-rate normalization intent.

**Caveat**: much of the surrounding kinetic-feedback / momentum-conserving-kick code path in `bh_update.c:270-296` is commented out — check whether `BH_FEEDBACK_ACTIVE` is actually enabled in current run configs before treating this as urgent.

---

## 4. Summary table

| Quantity | File:line | Status |
|---|---|---|
| Freefall time | `sfr_AGORA.c:160` | Missing `* cf_a3inv` on Density |
| Freefall time | `sfr_JEANS.c:198` | Missing `* cf_a3inv` on Density |
| Jeans length | `sfr_JEANS.c:62` | Missing `* cf_a3inv`; feeds SF criterion directly |
| Jeans mass | `sfr_JEANS.c:79` | Missing `* cf_a3inv`; feeds SF criterion directly |
| Freefall time | `individual_star_formation.c:94` | Missing `* cf_a3inv` on Density |
| Bondi relative velocity | `bh_accretion.c:351-353,489-490` | Missing `/ cf_atime` |
| Torque-accretion disk check | `bh_accretion.c:368-371,393` | Missing `/ cf_atime` (lower severity — sign test only) |
| BH thermal feedback energy | `bh_feedback.c:226`, `bh_update.c:287` | Missing `cf_atime^2` on Energy increment |
| Overdensity floor, HostHaloMass, halo-vel BH seeding, gravity/potential | various | Confirmed correct / not applicable |

---

## 5. Suggested next steps

1. Fix items (a) and (b) first — they're the same bug pattern repeated across three SF schemes, and (b) affects whether star formation triggers at all, not just its rate.
2. Confirm whether `BH_FEEDBACK_ACTIVE` / `TORQUE_ACCRETION` are enabled in current production `Config.sh` before prioritizing (c) and (d).
3. No files have been modified as part of this audit — all fixes are still open.

---

# Update 2026-10-23 (branch `grackle-batched-cooling`, working tree incl. uncommitted edits)

Static read-through only; nothing compiled or run.

## A. Status of the 2026-09-08 items
All of section 3 (a)-(d) are now fixed in the **uncommitted** working tree: freefall time (`sfr_AGORA.c`, `sfr_JEANS.c`, `individual_star_formation.c`), Jeans length/mass (`sfr_JEANS.c`), Bondi/torque relative velocity (`bh_accretion.c`, `/cf_atime`), BH thermal energy (`bh_update.c`, `* cf_atime^2`). Jeans-length fix re-derived: comoving L = c_s sqrt(pi a / (G rho_c)). Correct.

## B. New findings

### B1 (HIGH) Stellar age is a difference of scale factors, not of time
`src/stars/star_update.c:246`: `SP[i].Age = All.Time - SP[i].Birthtime;` then `star_age_yr = Age * cf_UnitTime_in_yr` (:251). With `ComovingIntegrationOn`, `All.Time` is `a`, so Age = Δa. `cf_UnitTime_in_yr` converts a *dln a* interval (it carries 1/H), so Age is neither cosmic time nor dln a. Δa = a·Δln a, so ages are too small by ~a (~7x at z=6) and set the wrong stellar lifetimes: delayed/mistimed winds, SNe and radiation, via `star_feedback_compute(dt, z, m, age)`.
(The timestep conversions at :105 and :250 are fine: dlna * cf_UnitTime = physical time.)
Fix: age_yr = `get_time_difference_in_Gyr(Birthtime, All.Time) * 1e9` in cosmological runs (helper exists, `timestep.c:702`; check its non-comoving branch), keep the existing path otherwise, and stop storing `Age` in mixed units.

### B2 (MEDIUM) SN terminal-momentum cap missing a factor of a
`src/stars/star_feedback.c:155`: `p_term /= (cf_UnitMass_in_Msun * cf_UnitVelocity_in_cm_per_s / 1e5)`. Code momenta are canonical (p = a p_phys), as the wind term (`star_interpolation.c:603`, `cf_UnitMomentum_in_cgs`, which has 1/a) already assumes. The cap is therefore too large by 1/a and rarely binds at high z, so the energy-conserving momentum from the quadratic is not limited as intended. Fix: divide by `cf_UnitMomentum_in_cgs` after converting Msun km/s to cgs.

### B3 (MEDIUM-LOW) SN lead time compared in the wrong time unit
`src/stars/star_density.c:550`: `sn_lead_time = SN_LeadTime / All.UnitTime_in_Megayears` (code *time*, no 1/H, no 1/h), but `TimeToSN` is converted with `/ cf_UnitTime_in_yr` in `star_interpolation.c:589` (dln a). Only matters with `TREE_BASED_TIMESTEPS && SUPERNOVAE`. Fix: `SN_LeadTime * 1e6 / All.cf_UnitTime_in_yr`.

### B4 (LOW, not expansion-factor) `h2` uninitialised in BH thermal feedback
`src/blackholes/bh_feedback.c:181` is commented out, but `if(r2 < h2)` at :213 still uses it unless `BH_CONSTANT_RADIUS` is defined. Also `hinv`/`u` are used at :214. Undefined behaviour in the default path. (`bh_jet_feedback.c` has the same `h2` use; not checked.)

### B5 (LOW / design) `BH_CONSTANT_RADIUS` is comoving
`All.BhRadius` is used directly as a comoving code length (`bh_accretion.c:331,360,565`, `bh_swallow.c:234,254`, `bh_feedback.c`), so the physical radius grows as a. Decide whether the parameter is meant to be physical; if so, divide by `cf_atime`.

### B6 (LOW, inert today) Jet feedback has no cosmology handling
`bh_jet_feedback.c:194` builds a physical energy and deposits it as `KineticFeed`; the kick/injection that would consume it is commented out in `bh_update.c`. If re-enabled it needs the a^2 energy factor and a canonical (a-scaled) kick.

### B7 (VERIFY) h factor in RT rates
`star_radiation.c:712,757,759,761`: `dt / cf_hubble_a / HubbleParam`. `cf_UnitTime_in_s` already carries 1/h, and Grackle `time_units = UnitTime/h`, so the extra `/HubbleParam` may double-count h (a factor-h rate error, not an `a` error). Worth confirming against the intended unit of the stored rates.

## C. Checked and correct
SN radius/energy conversions (`star_feedback.c:143-208`, via `cf_Unit*`); wind momentum (`star_interpolation.c:603`); radiation momentum (`star_radiation_tree.c:247-253`); star/BH energy injection (`star_update.c:443`: canonical energy); Grackle call (physical density `* cf_a3inv`, `a_value = cf_atime` each step, `comoving_coordinates = 0`); cooling `dtime = a*dt/(a*H)`; SF probability `dtime` (`starformation.c:193-202`).

## D. Suggested order
B1 first (silent, large), then B2/B3, then B4. B5-B7 are decisions or checks.

## E. Fix status (2026-10-23, later)
B1 (stellar age), B2 (SN momentum cap) and B3 (SN lead time) are fixed; see `documentation/source/cosmological_units.md` for the conventions. B4-B7 remain open. Nothing has been compiled or run in a comoving configuration yet.
