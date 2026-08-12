# Halo-Mass-Dependent Star Formation Threshold (`SF_THRESHOLD_HALO_MASS_DEPENDENT`)

## What it does

Raises the effective star-formation threshold for gas cells sitting in low-mass halos,
using `SphP[].HostHaloMass` (tagged on every gas cell at each on-the-fly FOF pass, see
`fof_halo_seeding.md`) rather than any new communication pass of its own.

This is a **single-knob, phenomenological stand-in**, not a first-principles calculation.
It's motivated by (but does not explicitly model) two different physical pictures that
both happen to point the same way — star formation should be *harder* in low-mass halos:

- a background UV/Lyman-Werner radiation field suppressing cooling/self-shielding in small
  halos (reionization-era feedback), or
- a Pop III -> Pop II transition, where the effective collapse threshold for metal-free gas
  differs from enriched gas.

Both are usually modelled with more machinery than this (a real radiation field, an actual
metallicity criterion) — this flag exists as a cheap, tunable proxy for either, using halo
mass as a single available handle, not as a replacement for either.

## Compile-time flag

`SF_THRESHOLD_HALO_MASS_DEPENDENT` in `Template-Config.sh` / `CosmoConfig.sh`, under
*Star Formation options*. Requires:

- `USE_SFR` and one of `EEOS_SF` / `AGORA_SF` / `JEANS_SF` (Makefile-enforced: `$(error
  SF_THRESHOLD_HALO_MASS_DEPENDENT requires USE_SFR...)`).
- `HALO_SEEDING` (C-level `#error` in `starformation.c`) — needs `SphP[].HostHaloMass`,
  which is only tagged when on-the-fly FOF is running.
- If `JEANS_SF` is the active scheme, also requires `JEANS_MASS_BASED` (C-level `#error`
  in `sfr_JEANS.c`) — the plain Jeans-length criterion has no threshold constant to
  modulate; see *Per-scheme application* below.

## Runtime parameters (`param.txt`)

Added under `#ifdef SF_THRESHOLD_HALO_MASS_DEPENDENT` in `src/io/parameters.c`:

- `MinHaloMassForNormalSF` — halo mass (code units) above which the ordinary threshold
  applies unchanged.
- `LowMassHaloThresholdFactor` — multiplier making star formation *harder* to trigger for
  gas cells whose host halo mass is below `MinHaloMassForNormalSF` (and above `0`).

Gas not currently in any FOF group (`HostHaloMass <= 0`) always uses the ordinary
threshold — this is deliberately an in-halo effect (self-shielding/collapse-threshold
physics inside a halo), not something meant to apply to general diffuse gas.

## Design: scheme-agnostic by construction

`EEOS_SF`, `AGORA_SF`, and `JEANS_SF` are mutually exclusive at compile time (only one may
be active — see the Makefile's `SF_MODELS` check) and each has its own, differently-named,
differently-unit'd threshold quantity (`PhysDensThresh`, `NumberDensThreshold`, a
Jeans-mass/length criterion). Rather than duplicating the halo-mass logic three times,
`src/star_formation/starformation.c` (always compiled under `USE_SFR`, regardless of which
scheme) provides one shared helper:

```c
double sf_threshold_halo_mass_factor(int i);
```

returning `All.LowMassHaloThresholdFactor` if cell `i`'s host halo mass is in
`(0, All.MinHaloMassForNormalSF)`, otherwise `1`. Each scheme's own star-formation
criterion calls this and folds it into its own threshold at the point of comparison — the
halo-mass dependence lives in exactly one place, and adding a fourth SF scheme in the
future only requires that scheme to call the same helper, not reimplement the logic.

## Related, independent gate: the comoving overdensity floor

`EEOS_SF` has always additionally required `All.ComovingIntegrationOn && dens >=
All.OverDensThresh` before treating a cell as star-forming, regardless of its own local
density/temperature gate — `AGORA_SF` and `JEANS_SF` did not, until this was added
alongside the halo-mass work above (same development pass, different motivation): with
less conservative local thresholds, both schemes could let diffuse,
cosmologically-unvirialized gas transiently satisfy their own criterion early in a run,
forming stars well before real structure had collapsed.

`All.OverDensThresh` is now set for all three schemes by one shared function,
`set_overdens_thresh()` in `starformation.c` (`All.OverDensThresh = All.CritOverDensity *
All.OmegaBaryon * 3 * All.Hubble^2 / (8*pi*G)`), called from `begrun.c` for
`EEOS_SF`/`AGORA_SF`/`JEANS_SF` alike. `CritOverDensity` (`param.txt`) is consequently now
required for all three schemes, not just `EEOS_SF` — the standard AREPO example value is
`57.7`. This is independent of, and unaffected by, `SF_THRESHOLD_HALO_MASS_DEPENDENT`
above — the two gates compose (both must pass), neither depends on the other being
enabled.

## Per-scheme application

- **`AGORA_SF`** (`sf_criteria()` in `sfr_AGORA.c`): `NumberDensThreshold` is multiplied by
  the factor before the `number_dens < threshold` gate. Straightforward — this scheme's
  inequality is oriented the same way the factor's semantics assume (factor `> 1` requires
  higher density, i.e. harder to form stars).

- **`EEOS_SF`** (`sfr_eEOS.c`): `PhysDensThresh` is multiplied by the factor, computed
  **per gas cell** (previously a single value computed once outside the per-particle loop
  in `cooling_and_starformation()` — moved inside the loop so it can vary per cell). Applied
  in both `cooling_and_starformation()`'s own gate and `get_starformation_rate()`'s
  independent copy of the same gate (the latter previously had a latent inconsistency where
  one comparison used a local `eos_dens_threshold` variable and a second, nearby comparison
  used `All.PhysDensThresh` directly — harmless while they were numerically identical, but
  fixed to use the same local variable now that they can differ).

  **Deliberately not modulated**: `calc_egyeff()`'s own internal use of `All.PhysDensThresh`
  (which sets the *shape* of the effective multiphase equation of state for a cell already
  identified as star-forming, via `rho = dmax(rho, PhysDensThresh)` and related ratios).
  Only the yes/no gate is halo-mass-dependent; once a cell is in the star-forming regime,
  its effective-EOS thermodynamics are unmodified. A collaborator wanting the whole eEOS
  curve to shift for low-mass-halo gas (not just the gate) would need to also pass a
  per-cell factor into `calc_egyeff()` (it already takes `i`, so `SphP[i].HostHaloMass` is
  available there too).

- **`JEANS_SF`** (`sf_criteria()` in `sfr_JEANS.c`): only meaningful with
  `JEANS_MASS_BASED`, where `JeansMassThreshold` is **divided** by the factor (not
  multiplied) before the `jeans_mass < threshold * cell_mass` gate — this scheme's
  inequality is the opposite orientation from the density-based schemes (star formation
  triggers when the Jeans mass is *smaller*), so dividing is what's needed to preserve the
  "factor `> 1` means harder to form stars" semantics across all three schemes. The plain
  Jeans-length criterion (`JEANS_SF` without `JEANS_MASS_BASED`) has no named threshold
  constant at all — just a hardcoded factor of `2.0` against the cell radius — so this
  combination is rejected at compile time rather than silently doing nothing.

## Known limitations / how to improve this

- **Hard step, not a smooth transition.** Star formation behaves identically for all halos
  above `MinHaloMassForNormalSF` and identically (but harder) for all halos below it, with
  a discontinuity at the boundary. A smooth interpolation (e.g. a sigmoid in
  `log(HostHaloMass)`) would avoid two similarly-massed halos straddling the cutoff having
  qualitatively different SF behaviour, at the cost of an extra transition-width parameter.
- **Halo mass only, not metallicity.** The Pop III/Pop II motivation is more standardly
  tied to a critical metallicity in the literature, not halo mass — this flag doesn't touch
  metallicity at all. A metallicity-aware version could combine `HostHaloMass` with each
  cell's own `SphP[].GasMetals` (already available, and already used by
  `BH_SEED_ON_ZERO_METALLICITY` elsewhere in the codebase) for a more physically-grounded
  criterion.
- **`EEOS_SF`'s effective-EOS shape is unmodified** for cells already past the gate (see
  above) — only whether a cell passes the gate at all is halo-mass-dependent.
- **Not run against a real halo catalogue.** Verified only by compiling successfully against
  all three SF schemes (including confirming the `JEANS_SF`-without-`JEANS_MASS_BASED`
  `#error` fires) — not tested in a live simulation. Recommend a small, cheap test run with
  a mix of halo masses straddling `MinHaloMassForNormalSF` before trusting it in production.
- **No feedback into the seeding channels documented in `fof_halo_seeding.md`** — this is a
  separate, independent use of `HostHaloMass`; seeding and star-formation-threshold
  behaviour don't interact with each other beyond sharing that one tagged field.
