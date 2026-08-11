# SOLAS Parameter Reference

Built directly from `src/io/parameters.c` (parameter names, required Config flags, and REAL/STRING/INT
type) and from grepping each parameter's usage site in `src/` for units and semantics, against `main`
as of Aug 2026. This supersedes the parameter tables in `Parameters.md`, `Code.md`, and the
Development Guide's §5 flag table — please point people here and correct the others to link in
rather than re-listing values, so there's one place this drifts out of sync from instead of four.

**Legend:** ✅ verified against a source usage site · ⚠️ inferred from context, not directly confirmed
· ❓ genuinely unknown from source alone — needs input from whoever wrote/last touched it.

**Known gap:** no committed example run turns on `AGORA_SF`, `JEANS_SF`, `BLACKHOLES`, `RADIATION`,
or `HALO_SEEDING` (both shipped examples use plain `EEOS_SF`), so there is no in-repo "known-good"
value for most of the parameters below. Where a real production value exists, it's quoted with its
source. Everywhere else, "Recommended" is blank rather than invented — please fill these in as the
team settles on values, ideally by adding a worked example config rather than only editing this file.

**Where this should actually live:** `documentation/source/` already has a working Sphinx setup —
`index.md`'s `toctree` wires in `parameterfile.md`/`config-options.md` (the untouched public-AREPO
docs, which is why none of this is in there yet), plus `technical-reference.md` and
`fof_halo_seeding.md` (SOLAS-specific, already merged into `main`, already in the toctree). This
file fits that pattern as a new `documentation/source/solas-parameters.md`, added to the `toctree`.
One more thing while you're in there: `sf_threshold_halo_mass.md` exists in `documentation/source/`
but is **not** in `index.md`'s toctree — it's built but orphaned, unreachable from the docs
homepage. Worth a one-line fix alongside this.

---

## Star formation — general

| Parameter | Required flag(s) | Type | Units | Recommended | Notes |
|---|---|---|---|---|---|
| `CritOverDensity` | `EEOS_SF` \| `AGORA_SF` \| `JEANS_SF` (needs `USE_SFR`) | REAL | ✅ dimensionless — `documentation/source/parameterfile.md` gives the exact formula: `rho_th = CritOverDensity × 3 × Omega_b × H² / (8πG)` (redshift-independent) | `57.7` (both shipped `EEOS_SF` examples) | Applies across all three SF flags per the `#if` guard in `parameters.c` — not Jeans-specific. **Used for comoving (cosmological) runs**; `CritPhysDensity` is the non-comoving equivalent below — they're alternates for the same role, not both active at once. |
| `CritPhysDensity` | same as above | REAL | ✅ cm⁻³, physical — **used instead of `CritOverDensity` for non-comoving runs**, per `documentation/source/parameterfile.md` | `0` in both shipped examples — including the non-comoving `galaxy_merger` one, which per the doc should be relying on this one, not `CritOverDensity`. Worth a quick check with Chris/Balu on whether that example's `param.txt` is actually doing what it intends. | |

## `EEOS_SF` (public AREPO's default Springel & Hernquist 2003 model — not a SOLAS addition, included for completeness since it's what both shipped examples actually use)

| Parameter | Type | Units | Recommended | Notes |
|---|---|---|---|---|
| `MaxSfrTimescale` | REAL | ✅ internal time units (per in-repo comment) | `2.27` | |
| `TempSupernova` | REAL | ✅ Kelvin (per in-repo comment) | `5.73e7` | |
| `TempClouds` | REAL | ✅ Kelvin (per in-repo comment) | `1000.0` | |
| `FactorEVP` | REAL | ⚠️ dimensionless | `573.0` | |
| `FactorSN` | REAL | ⚠️ dimensionless (mass fraction) | `0.1` | |
| `TemperatureThresh` | REAL | ⚠️ Kelvin, by analogy with `TempSupernova`/`TempClouds` | `1e6` | |

## `AGORA_SF`

| Parameter | Type | Units | Notes |
|---|---|---|---|
| `NumberDensThreshold` | REAL | ✅ **cm⁻³, physical** — `sfr_AGORA.c` computes `number_dens = SphP[i].Density * cf_UnitDensity_in_cgs / mu / PROTONMASS` and compares directly against this | Same name and role as `INDIVIDUAL_STAR_BY_STAR_FORMATION`'s threshold below — they're the same parameter tag, read once, shared by whichever flag is active. |
| `TemperatureThreshold` | REAL | ✅ Kelvin — compared directly against a temperature computed in K | SF disabled above this temperature. |
| `StarFormationEfficiency` | REAL | ✅ dimensionless — `SFR = efficiency × Mass / t_freefall` | Same tag reused by `JEANS_SF` and `INDIVIDUAL_STAR_BY_STAR_FORMATION`. |

## `SF_THRESHOLD_HALO_MASS_DEPENDENT` (found in source — not mentioned on any current wiki page)

Multiplies `NumberDensThreshold` by a halo-mass-dependent factor (`sf_threshold_halo_mass_factor()`)
— this is presumably the mechanism behind the `sf_threshold_halo_mass.md` doc referenced from
`Home.md`'s FOF-seeding section, but that doc lives on the (now-merged) `merge-fof-into-star-feedback`
branch and I haven't read it. Worth linking once someone confirms the current path.

| Parameter | Type | Units | Notes |
|---|---|---|---|
| `MinHaloMassForNormalSF` | REAL | ❓ presumably code mass units, unconfirmed | |
| `LowMassHaloThresholdFactor` | REAL | ❓ presumably dimensionless multiplier | |

## `JEANS_SF`

| Parameter | Required flag(s) | Type | Units | Notes |
|---|---|---|---|---|
| `JeansMassThreshold` | `JEANS_SF` + `JEANS_MASS_BASED` only | REAL | ✅ dimensionless — SF enabled if `jeans_mass < JeansMassThreshold × cell_mass` | Not read at all in plain `JEANS_SF` (length-based) mode. |
| `StarFormationEfficiency` | `JEANS_SF` (either mode) | REAL | ✅ dimensionless, same as `AGORA_SF` above | |

## `INDIVIDUAL_STAR_BY_STAR_FORMATION`

Reads the same three tags as `AGORA_SF` (`NumberDensThreshold`, `TemperatureThreshold`,
`StarFormationEfficiency`) — same units, same meaning, just gating a different formation mode
downstream. Additionally:

| Parameter | Required flag(s) | Type | Units | Notes |
|---|---|---|---|---|
| `IMF` | `STAR_PARTICLES` | **INT**, not string | ❓ presumably an IMF-selector enum | The old `Parameters.md` implied this was a path/name; source says `INT`. `technical-reference.md` mentions stars are "IMF-binned across `NBINS = 114` mass bins" but doesn't give the integer→IMF-model mapping — confirm with whoever wrote `star_particle.c`. |
| `StarTablesFile` | `STAR_FEEDBACK_ACTIVE` | STRING | path | Stellar-property lookup table, as in the Development Guide. |
| `SN_LeadTime` | `TREE_BASED_TIMESTEPS` + `SUPERNOVAE` | REAL | ❓ presumably code time units | Only read when both flags are set together — easy to silently omit. |

## Stellar feedback

| Parameter | Required flag(s) | Type | Units | Notes |
|---|---|---|---|---|
| `SN_HostShellSweepFrac` | `SUPERNOVAE` | REAL | ✅ dimensionless fraction, **internally capped at 0.9** regardless of input (`fmin(All.SN_HostShellSweepFrac, 0.9)` in `star_feedback.c`) | Note the underscore — `SN_HostShellSweepFrac`, not `SNHostShellSweepFrac` as I'd written earlier. |

## Radiation

| Parameter | Required flag(s) | Type | Units | Notes |
|---|---|---|---|---|
| `RaySplitFactor` | `STAR_RADIATION_ACTIVE` | REAL | ✅ dimensionless multiplier on node child count (`star_radiation_tree.c`) | |
| `RadOpeningAngle` | `RAD_OPENING_ANGLE` | REAL | ⚠️ dimensionless angle, compared as `(A_proj/dist²) < RadOpeningAngle²` — same style as AREPO's gravity opening criterion, so likely radians by convention but not explicitly confirmed | |
| `NodeAspectRatio` | `RAD_OPENING_ANGLE` | REAL | ❓ | |
| `RTIonizationTimestepFraction` | `PHOTOIONIZATION` | REAL | ✅ dimensionless fraction (`eps_ion` in `star_radiation.c`) | |
| `IRDtauMomentumBoostCoeff` | `RADIATION_PRESSURE` | REAL | — | **Dead parameter as of this checkout** — read into `All.IRDtauMomentumBoostCoeff` but never referenced anywhere else in `src/`. Flag to the team before anyone tunes it. |

## Black holes — density/radius

| Parameter | Required flag(s) | Type | Units | Notes |
|---|---|---|---|---|
| `BhRadius` | `BH_ACTIVE` + `BH_CONSTANT_RADIUS` | REAL | ✅ code length units — used directly in `r² < BhRadius²` distance comparisons | |
| `BhDesNgb` | `BH_ACTIVE`, else-branch (adaptive radius) | REAL | ✅ dimensionless, mass-weighted neighbour count (compared against `NgbsMass / TargetGasMass`) | Distinct from gas's `DesNumNgb` and BH-jet's separate check in `bh_jet_density.c`. |
| `BhDesDev` | same as above | REAL | ✅ dimensionless, same convention | |
| `HMaxFactor` | same as above | REAL | ✅ dimensionless — multiplies `SofteningTable[...]` to get a max search radius | |

## Black hole accretion

| Parameter | Required flag(s) | Type | Units | Notes |
|---|---|---|---|---|
| `Epsilon_r` | `BH_ACCRETION_ACTIVE` | REAL | ⚠️ dimensionless radiative efficiency (standard accretion-physics usage, e.g. `energy = Epsilon_f × Epsilon_r × Ṁc²`) | |
| `Epsilon_T` | `TORQUE_ACCRETION` | REAL | ❓ normalisation coefficient for the torque accretion rate — units depend on the rest of the formula in `bh_accretion.c:618`, not fully traced here | |
| `ADP_tvisc` | `ADP_ACCRETION` | REAL | ⚠️ code time units — `mdot_visc = M_disc / ADP_tvisc` | |
| `ADP_tcap` | `ADP_ACCRETION` | REAL | ⚠️ code time units — `dM_to_disc = M_res × (bh_timestep / ADP_tcap)` | |
| `ADP_EddFactor` | `ADP_ACCRETION` | REAL | ✅ dimensionless — `mdot_cap = ADP_EddFactor × EddingtonRate` | |

## Black hole feedback

| Parameter | Required flag(s) | Type | Units | Notes |
|---|---|---|---|---|
| `Epsilon_f` | `BH_FEEDBACK_ACTIVE` | REAL | ⚠️ dimensionless coupling efficiency | Also used by `bh_jet_feedback.c`, so applies to jets too even though it's only gated on `BH_FEEDBACK_ACTIVE`. |
| `Mload` | `BH_FEEDBACK_ACTIVE` | REAL | ⚠️ dimensionless, expected in [0,1] — appears as both `Mload` and `(1 − Mload)` in the mass/energy split | |

## Black hole refinement

`REFINEMENT_AROUND_BH` has **three** distinct sub-modes reading different parameter sets — the
Development Guide's "==0 fixed radius, ==1 hybrid" description undersells this; there's a third,
default branch:

| Sub-mode | Parameters | Notes |
|---|---|---|
| `REFINEMENT_AROUND_BH_FIXED` | `RefBHRadius`, `RefBHMinCellRadius`, `RefBHMaxCellRadius` | |
| `REFINEMENT_AROUND_BH_HYBRID` | `RefBHRadiusHSML`, `RefBHMinCellRadius`, `RefBHMaxCellRadius` | |
| default (neither defined) | `RefBHRadiusHSML`, `RefBHMinCellRadiusRBondi`, `RefBHMaxCellRadiusHSML` | Note the parameter **names differ** from the hybrid mode despite `RefBHRadiusHSML` being shared — `RefBHMinCellRadius` vs `RefBHMinCellRadiusRBondi` are not interchangeable. |

Always read under `REFINEMENT_AROUND_BH` regardless of sub-mode:

| Parameter | Type | Units | Notes |
|---|---|---|---|
| `RefBHMinCellMass` | REAL | ❓ presumably code mass units | |
| `RefBHLowerFactorC` | REAL | ❓ | |

Separately, star-hosting-cell refinement is gated by **`STAR_HOST_REFINEMENT`** — not bundled
under `REFINEMENT_AROUND_BH` and not currently in the Development Guide's flag table at all:

| Parameter | Required flag | Type | Units |
|---|---|---|---|
| `RefStarsPerCell` | `STAR_HOST_REFINEMENT` | REAL | ❓ presumably a count threshold |

## Black hole merging (`BH_MERGER` — not mentioned in the Development Guide or Code page at all)

Full mechanism (donor/survivor selection, drain clamp, merge-count bookkeeping) is documented in
`documentation/source/technical-reference.md` §2.3 (`bh_merger.c`) — worth reading before touching
either parameter below rather than inferring behaviour from the names.

| Parameter | Type | Units | Notes |
|---|---|---|---|
| `BhMergerRadiusFactor` | REAL | ⚠️ dimensionless — `merger_radius = BhMergerRadiusFactor × h` | |
| `BhMergerRadiusCriterion` | STRING | — | e.g. `"HSML"` reproduces old default behaviour exactly, per the comment at `parameters.c:760` — see `enum bh_merger_radius_criterion` in `src/blackholes/bh.h` for the other options. **This is a mandatory tag under `BH_MERGER`** — existing param files need it added explicitly. |

## On-the-fly FOF halo / BH seeding (`HALO_SEEDING`)

| Parameter | Required flag(s) | Type | Units | Notes |
|---|---|---|---|---|
| `TimeOfFirstHaloFinding` | `HALO_SEEDING` (requires `FOF`) | REAL | ✅ code time units (same convention as `All.Time`; scale factor for cosmological runs) | |
| `TimeBetweenHaloFinding` | `HALO_SEEDING` | REAL | ✅ dimensionless multiplicative factor (`NextTimeOfHaloFinding *= this`) | Geometric, not additive — a value like `1.1` compounds each check. |
| `MinHaloMassForFOFSeeding` | `BH_SEED_ON_MASS` | REAL | ✅ code mass units — compared directly to FOF `Group[n].Mass` | |
| `ZeroMetallicityThresholdForFOFSeeding` | `BH_SEED_ON_ZERO_METALLICITY` | REAL | ✅ dimensionless metal mass fraction — compared to `Group[n].MaxGasMetallicity` | |
| `MinVelDispForFOFSeeding` | `BH_SEED_ON_VELDISP` | REAL | ✅ km/s-equivalent DM 3D velocity dispersion (`sigma_3D`, mass-weighted, from `Group[].VelDispDM`) | Not on the **wiki's** `FOF-seeding-branch` page, but fully documented in `documentation/source/fof_halo_seeding.md` (already in `main`) — link to that instead of duplicating; the wiki page just needs a pointer added, same pattern as `Technical-Reference.md`. |
| `PotentialDonorSearchNSoft` | `BH_SEED_ON_POTENTIAL_POSITION` (also requires `EVALPOTENTIAL`, not just `BLACKHOLE_SEEDING` — correction to my earlier note) | REAL | ✅ multiple of the DM gravitational softening length | **Suggested default 3–5**, explicitly flagged in `fof_halo_seeding.md` as "untuned against real halos, worth revisiting" — don't treat that range as validated. |
| `BlackHoleSeedMass` | `BLACKHOLE_SEEDING` | REAL | ✅ code mass units — compared directly to `P[igas].Mass` | |

## Cooling / metals

| Parameter | Required flag | Type | Units | Notes |
|---|---|---|---|---|
| `InitMetallicityinSolar` | `METALS` | REAL | ⚠️ solar units (per name) | |
| `TreecoolFile` | `COOLING` | STRING | path | |
| `GrackleDataFile` | `USE_GRACKLE` | STRING | path | |

---

## How this was built, so it can be regenerated

```bash
git clone https://github.com/solas-sims/arepo_solas.git
cd arepo_solas
grep -n 'strcpy(tag\[nt\]' src/io/parameters.c        # full parameter list + line numbers
grep -n 'All\.<ParamName>\b' src/**/*.c                # usage sites, for units/semantics
```
Everything marked ✅ above was confirmed this way. Everything marked ❓ needs five minutes with the
same grep against whoever's local checkout has the relevant physics turned on and running — I'd
rather leave a gap flagged than fill it with a plausible-sounding guess.
