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
a true no-op); all 19 standalone rate-law tests pass.

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
3. **Rate-law constants are literature-standard placeholders, not verified
   against a specific source.** Every constant in `src/dust/dust.h` is
   commented `VERIFY` — none have been checked against Li, Narayanan & Dave
   (2019) directly.

## What needs to be done

### Before any physics results can be trusted
- [ ] **Verify rate-law constants** in `src/dust/dust.h`
  (`DUST_SN_CONDENSATION_EFFICIENCY`, `DUST_GROWTH_TAU_REF_GYR`/`REF_NH_CGS`/`REF_TEMP_K`,
  `DUST_SPUTTER_TAU_REF_GYR`/`REF_NH_CGS`/`REF_TEMP_K`/`OMEGA`) against
  Li, Narayanan & Dave (2019) directly, with collaborators.
- [ ] **Decide on the wind-channel dust production gap.** Either accept SN-only
  as a standing Phase 1 limitation, or prioritize splitting the `WINDS`
  yield table into AGB vs. massive-star components (a `src/stars/`
  change, not a dust-module change).

### To actually run and validate Phase 1
- [ ] Supply an isolated-disk (or similar idealized ISM) **initial conditions
  file** for `examples/dust_isolated_disk_3d/`.
- [ ] Supply a **Grackle chemistry/cooling data file** (`GrackleDataFile` in
  `param.txt`).
- [ ] Supply the **stellar yield/feedback table** (`StarTablesFile` in
  `param.txt`, consumed by `load_star_tables()`).
- [ ] Tune the placeholder parameters in `param.txt` (`BoxSize`,
  `ReferenceGasPartMass`, `MeanVolume`, softening lengths, `SN_LeadTime`,
  `SN_HostShellSweepFrac`, `NumberDensThreshold`, `TemperatureThreshold`,
  `StarFormationEfficiency`, `InitMetallicityinSolar`) to match whatever IC
  is supplied.
- [ ] Run the **Phase 1 exit gate**:
  - Zero-rate bit-for-bit check: build with `DUST` on and off, with rate
    constants effectively zeroed (e.g. `DUST_SN_CONDENSATION_EFFICIENCY 0`),
    confirm the two runs agree to machine precision outside `DustMass`.
  - Metal mass conservation: `GasMetals` (gas-phase + dust-phase) conserved
    over the run.
  - D/Z–metallicity trend qualitatively consistent with Li et al. (2019)'s
    published relation.
- [ ] Note this build environment's local `Config_AGORA.sh`/`USE_GRACKLE`
  path currently fails to compile against the Grackle library installed on
  this dev machine (API version mismatch, pre-existing and unrelated to
  dust) — verification here used a `USE_GRACKLE`-free config instead. Confirm
  the real run environment (cluster) has a matching Grackle version before
  using `Config_AGORA.sh`-style flags for real.

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
