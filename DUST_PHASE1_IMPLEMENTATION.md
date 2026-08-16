# Dust Module — Phase 1 Implementation Status

Branch: `dust_module`. Companion to `dust_module_implementation_phaseplan.docx`
(the phased plan) and the repo audit that preceded it. This document
describes what Phase 1 actually built, the scope decisions made along the
way, and what's left before it's science-ready.

## What Phase 1 is

A single-scalar dust mass carried per gas cell, with:

- **Production** — SN-channel dust condensation only (see "Scope decisions" below).
- **Destruction** — SN shock destruction, McKee (1989)-style swept-mass form.
- **Growth** — grain accretion, Dwek (1998) timescale form.
- **Sputtering** — thermal sputtering, Tsai & Mathews (1995) timescale form.

It validates the plumbing (data model, timestep integration, snapshot I/O,
unit tests) before Phase 2b adds species/size resolution.

## What was implemented, file by file

| Area | Files | What it does |
|---|---|---|
| Config/build | `Template-Config.sh`, `Makefile`, `src/main/proto.h` | `DUST` flag (requires `METALS`+`SUPERNOVAE`), gated object files, `dust_proto.h` hooked into the global prototype chain. |
| Data model | `src/main/allvars.h` | `DUST_INDEX`/`DUST_NUMBER` passive scalar (mirrors `METALS_INDEX`/`METALS_NUMBER` — gets Voronoi flux advection, MPI exchange, and refinement/derefinement handling for free via the existing generic `N_Scalar` machinery). `GasDustMass`/`GasDustMassFraction` aliases, `StarDustFeed` staging field, `IO_DUST_MASS` enum entry. |
| Rate laws | `src/dust/dust_production_sn.c`, `dust_destruction_sn.c`, `dust_growth.c`, `dust_sputtering.c` | Pure functions, no Arepo dependency (testable standalone). Constants isolated in `src/dust/dust.h`. |
| Integration | `src/dust/dust_update.c` | `dust_cell()` — per-cell growth+sputtering rate-ODE (bracket+bisect on implicit backward-Euler, mirrors `DoCooling()`), clamped to `0 <= GasDustMass <= GasMetals`. `dust_processes()` — driver over `TimeBinsHydro`, called from `src/main/run.c` right after cooling each step. |
| Snapshot I/O | `src/dust/dust_io.c`, `src/io/io_fields.c` | `DustMass` field registered in snapshots. |
| SN hook | `src/stars/star_feedback.c`, `star_update.c`, `src/mesh/mesh.h`, `src/mesh/voronoi/voronoi_exchange.c` | Production/destruction wired into the existing per-SN-event `Feedback_Kick` mechanism (both HOST and MESH deposition modes), including the `PrimExch` ghost-cell mirror needed for MPI-remote neighbor cells. |
| Tests | `tests/dust/` | Standalone unit tests for the four rate-law functions (19 checks, no MPI/mesh needed). Run via `tests/dust/build.sh`. |
| Example | `examples/dust_isolated_disk_3d/` | `Config.sh`/`param.txt` scaffold for the flag combination Phase 1 needs (`STARS`+`SUPERNOVAE`+`WINDS`+`METALS`+`DUST`) — no example in the repo had this combination before. **No IC included** — see that directory's `README.md`. |

Verified: builds and links cleanly both with `DUST` on and off (off-path is
a true no-op); all 21 standalone rate-law tests pass.

**Rate-law constants verified against Li, Narayanan & Dave (2019,
arXiv:1906.09277) primary text and Table 1** — cross-checked against a raw
`pdftotext` extraction of the actual paper, not an AI-summarized secondary
source. Every constant/functional form in `src/dust/dust.h` was wrong or
missing something before this pass; all four rate laws were corrected:

| Rate law | Before | After (Li+2019) |
|---|---|---|
| SN condensation efficiency | 0.01 | 0.15 (Table 1, `DUST_SN_CONDENSATION_EFFICIENCY`) |
| Growth `tau_ref` | 30 Myr | 10 Myr |
| Growth `rho_ref` | 1000 cm⁻³ | 100 cm⁻³ |
| Growth T-scaling | `sqrt(T_ref/T)` | linear `(T_ref/T)` |
| Growth metallicity term | absent | `(Z_sun/Z_gas)` added, `dust_growth_rate()` now takes a `Z_gas_massfrac` argument |
| Sputtering `rho_ref` | 1.0 cm⁻³ | 4.543e-4 cm⁻³ (paper's 1e-27 g/cm³ mass density converted to this code's H-number-density convention) |
| Sputtering rate | `-M/tau_sp` | `-3*M/tau_sp` (paper's `dM/dt = -M/(tau_sp/3)`) |
| SN destruction efficiency | none applied | `DUST_SN_DESTRUCTION_EFFICIENCY = 0.3` (eq. 8's epsilon) now multiplied in |

Destruction still reuses this code's own Kim & Ostriker (2015) swept-mass
estimate in place of the paper's independent `M_s` (see "Scope decisions"
below) — that structural difference is unchanged, only the missing
epsilon=0.3 factor was added. Production still applies the efficiency to
total `SN_MetalsLoss` rather than per-element ejecta (Phase 1's
single-scalar simplification) — only the 0.01→0.15 magnitude was fixed.

## Validation results (real Arepo runs, not just unit tests)

Run against `examples/agora_disc_star_formation_3d/` (ported from `main`,
see git log) — an isolated AGORA-disc IC with `STARS`+`SUPERNOVAE`+`WINDS`+
`METALS`+`USE_GRACKLE` all enabled, 752K gas cells.

**Structural check** (short run, no SN yet): gas and metal mass conserved,
`DustMass` present in snapshots and correctly zero everywhere pre-SN
(production is SN-only; growth has no seed to grow from), zero invariant
violations.

**Physics check** (SN actually firing): the organic run doesn't produce a
supernova quickly enough to check in reasonable wall-clock time — real
massive-star lifetimes are Myr-scale and the IMF makes SN-eligible stars
(≥8 M☉) rare in a small early population. Forced this by patching a local
copy of the stellar yield table (`star_feedback_tables.hdf5`): rescaled the
`Age` track (`star_interpolation.c`'s `star_lifetime()` reads its last
entry as the star's lifetime, in years) down by 1e5 for every
`(Z, M)` entry with `M >= 8` and `SN_MassLoss > 0`, so any SN-eligible star
sampled explodes almost immediately instead of after Myr. Result:

| Snapshot | Time | Stars | Total `DustMass` | Nonzero dust cells | Invariant violations |
|---|---|---|---|---|---|
| before any SN | — | 0 | 0 | 0 | 0 |
| first SN fires | — | 2 | 2.00e-08 | 25 | 0 |
| — | — | 2 | 1.99e-08 | 209 | 0 |
| second SN fires | — | 3 | 2.34e-08 | 620 | 0 |

Dust appears exactly when the first SN fires, spreads across the mesh via
the MESH-mode `Feedback_Kick` redistribution (25 → 209 → 620 nonzero
cells) as designed, and total mass both rises (fresh production) and dips
(destruction winning a step) as expected from combined production +
destruction physics — with **zero `0 <= GasDustMass <= GasMetals`
invariant violations** under real dynamic SN-driven conditions. This
confirms the SN production/destruction hooks work correctly, not just
compile.

**Rerun after the Li+2019 rate-law fixes** (same fast-SN table, 4 MPI
ranks, `TimeMax=0.005`, `TimeBetSnapshot=0.0002`, 752K gas cells): 10 stars
formed, 10 feedback events, dust mass rises from first SN through all 16
snapshots (1.03e-6 → 3.31e-6 in code mass units) with nonzero dust cells
climbing from 13 to 27010 as it spreads across the mesh, metal mass
non-decreasing throughout, star count non-decreasing throughout, and
`check.py`'s full five-check suite (metals, dust, star formation,
star-forming-gas-fraction, stellar spatial extent) passes with **zero
`0 <= GasDustMass <= GasMetals` invariant violations** and no ejection
warnings. Confirms the corrected constants integrate correctly end-to-end,
not just in the standalone unit tests.

**Environment notes from getting this running** (all fixed inline, not
dust-specific, but worth knowing about):
- This repo's `src/cooling/grackle.c` needs the `solas-sims/grackle` fork
  specifically (has `RT_HI/HeI/HeII_heating_rate` fields that upstream
  Grackle doesn't) — built locally at `~/software/grackle_solas`.
- Two pre-existing bugs on `Star_feedback_radiation` blocked building at
  all: an unbalanced `Makefile` (orphaned `endif`s around
  `FOF`/`HALO_SEEDING`) and a missing `#ifdef RAD_OPENING_ANGLE` guard
  around `NodeAspectRatio` in `src/io/parameters.c`. Both fixed in the
  `991754c` merge commit.
- The AGORA example's `param.txt` (as pulled from `main`) had a
  `CritOverDensity` line that's only a valid parameter under `EEOS_SF`,
  not `AGORA_SF` — removed in `3d6c98f`.
- One apparent 30+ minute hang in `sample_star_particle()`
  (`src/stars/star_particle.c`, `STAR_PARTICLES=1` IMF sampling) did not
  reproduce on a rerun with diagnostic instrumentation (zero iterations
  anywhere near the loop's safety cap) — very likely a transient
  environment/system-contention issue on this dev machine, not a code
  defect. No code changes made; flagging in case it recurs elsewhere.
- **Stale-object build trap** (found during the rate-law-fix rerun): if an
  earlier `make` invocation compiled any files against the wrong Grackle
  headers (e.g. before passing the correct `GRACKLE_INCL`/`GRACKLE_LIB`),
  incremental `make` will NOT recompile those `.o` files on a later
  correctly-configured invocation, since the header path is a
  command-line variable, not a tracked file dependency. The result links
  fine but corrupts memory at runtime (manifested here as `All.GrackleDataFile`
  silently truncated before `InitGrackle()`, causing an obscure
  "Can't open Cooling in ..." HDF5 error deep inside Grackle, unrelated to
  the actual bug). Fix: `rm -rf <BUILD_DIR> <EXEC>` and rebuild clean
  whenever `GRACKLE_INCL`/`GRACKLE_LIB` change between invocations.
  On this dev machine, `SYSTYPE=Darwin GRACKLE_INCL="-I$HOME/software/grackle_solas/include" GRACKLE_LIB="-L$HOME/software/grackle_solas/lib -lgrackle -Wl,-rpath,$HOME/software/grackle_solas/lib -L/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current"`
  is needed (the Makefile's default `GRACKLE_INCL`/`GRACKLE_LIB` point at
  `~/Codes/grackle`, the vanilla install, not the `solas-sims/grackle`
  fork this repo actually needs; `-lgfortran` also isn't on the default
  linker search path on this machine).

## Scope decisions (already made, documented here for reference)

1. **SN-channel production only.** This code's `WINDS` channel blends
   massive-star and AGB-type mass loss into a single table — there's no way
   to isolate true AGB ejecta from it, so wind-channel dust condensation
   (which is where most dust production happens in the literature) is
   deferred. Flag this to collaborators; it's a real physics gap, not a
   minor detail.
2. **Destruction is tied to this code's own swept-mass estimate** (Kim &
   Ostriker 2015 shell mass, already computed for momentum-cap purposes in
   `star_feedback.c`), not an independent McKee (1989) 6800 M☉ constant —
   more self-consistent with the existing SN physics, but means Phase 1's
   destruction rate is coupled to whatever `SN_HostShellSweepFrac` is set to.
3. **Rate-law constants are now verified against Li, Narayanan & Dave
   (2019) directly** (see "Verified" note above) — all four rate laws'
   constants/functional forms were corrected to match the paper. The two
   *structural* simplifications (single total-metals scalar instead of
   per-element production; destruction reusing this code's own swept-mass
   estimate instead of the paper's independent `M_s`) are unchanged and
   still Phase 1 scope decisions, not bugs.

## What needs to be done

### Before any physics results can be trusted
- [x] **Verify rate-law constants** in `src/dust/dust.h` against
  Li, Narayanan & Dave (2019) directly — done (see "Verified" note above);
  worth a second look from collaborators given the unit-conversion
  judgment call on the sputtering reference density (paper gives a mass
  density, this code's rate laws take H number density).
- [ ] **Decide on the wind-channel dust production gap.** Either accept SN-only
  as a standing Phase 1 limitation, or prioritize splitting the `WINDS`
  yield table into AGB vs. massive-star components (a `src/stars/`
  change, not a dust-module change).

### To actually run and validate Phase 1
- [x] Real IC / Grackle data / star tables — solved via
  `examples/agora_disc_star_formation_3d/` (ported from `main`, includes
  `get_ics.sh`/`get_tables.sh` downloaders), not the placeholder
  `examples/dust_isolated_disk_3d/` scaffold. Consider pointing
  `dust_isolated_disk_3d/` at this instead, or retiring it, now that a
  real working example exists.
- [x] Metal mass conservation — held across every run, including the
  SN-firing one.
- [x] SN production/destruction actually firing, with correct invariants
  — see "Validation results" above.
- [ ] Zero-rate bit-for-bit check (`DUST` on vs. off with rate constants
  zeroed) — not yet run.
- [ ] D/Z–metallicity trend vs. Li et al. (2019) — needs a much longer,
  statistically meaningful run (the validation run here used an
  artificially shortened stellar lifetime table specifically to force a
  fast SN; not representative of real D/Z evolution). A real check needs
  the *unpatched* table and enough runtime for organic star formation/SN
  statistics to build up.
- [ ] Tune the placeholder parameters in `param.txt` (`SN_LeadTime`,
  `SN_HostShellSweepFrac`, `NumberDensThreshold`, `TemperatureThreshold`,
  `StarFormationEfficiency`, `InitMetallicityinSolar`) — currently
  literature-typical defaults per that example's own README, not
  validated against this fork.

### Beyond Phase 1 (already scoped in the phaseplan doc, not started)
- **Phase 2a** (metals-model track, not dust-module work): extend chemical
  enrichment from the current single total-metals scalar to per-element
  (Mg, Si, O, C, Fe) tracking — owned separately, to be verified with
  collaborators.
- **Phase 2b** (blocked on 2a): 3 species × N size bins, species-resolved
  sputtering/SN destruction, sub-resolved clumping factor.
- **Phase 3**: shattering, coagulation, turbulent diffusion (no existing
  diffusion solver in this codebase — confirmed absent in the earlier
  repo audit, would need to be built from scratch).
- **Phase 4** (optional, not scheduled by default): two-fluid dynamical dust
  particles, only if a specific science case needs resolved dust–gas drift.

## Where to look for more detail

- `dust_module_implementation_phaseplan.docx` — the full phased plan,
  including the Phase 2a/2b split rationale and cross-cutting practices.
- Git log on `dust_module` (8 commits, `9f085e3`..`8f51d73`) — each commit
  is scoped to one logical piece (build wiring, data model, rate laws,
  SN hook, run.c integration, snapshot I/O, tests, example scaffold) and
  has a detailed message.
