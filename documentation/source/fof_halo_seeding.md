# On-the-fly FOF Halo Seeding (`fof` branch)

This document describes the work merged from the `arepo_solas_features` repository's `fof`
branch into `Star_feedback_radiation` (merge commit, July 2026). It covers 29 commits of
development and adds a mechanism for seeding halos — currently used for black hole seeding —
during the code's regular on-the-fly friends-of-friends (FOF) passes.

## What it does

Every time the on-the-fly FOF group finder runs, each identified halo is now optionally
checked against seeding criteria. If a halo qualifies, a seed event is recorded; seed events
are later converted into new black hole particles at the next safe synchronisation point in
the timestep loop (mirroring how star particles are spawned).

This replaces an older, simpler mechanism (`FIND_HALOS`) that lived on `Star_feedback_radiation`
and triggered `fof_seeding()` on a fixed wall/simulation-time schedule, independent of the FOF
pass. That mechanism has been retired as part of the merge — see "Note on the merge" below.

## Compile-time flags

Set in `Template-Config.sh` / `CosmoConfig.sh`, under *Inline halo finding*:

| Flag | Requires | Purpose |
|---|---|---|
| `HALO_SEEDING` | `FOF` | Enables tracking of halo properties needed for seeding (densest gas cell, host halo mass, optionally metallicity/velocity dispersion/potential minimum) during each FOF pass, and the seed-event list mechanism. |
| `BLACKHOLE_SEEDING` | `HALO_SEEDING` | Enables spawning black hole particles from seed events. Requires at least one seeding channel below. |
| `BH_SEED_ON_MASS` | `BLACKHOLE_SEEDING` | Seeding channel: seed once a halo's FOF mass exceeds `MinHaloMassForFOFSeeding`. |
| `BH_SEED_ON_ZERO_METALLICITY` | `BLACKHOLE_SEEDING`, `METALS` | Seeding channel: seed "pristine" halos, i.e. halos whose most metal-enriched gas cell still has a metal mass fraction at or below `ZeroMetallicityThresholdForFOFSeeding`. |
| `BH_SEED_ON_VELDISP` | `BLACKHOLE_SEEDING` | Seeding channel: seed once a halo's DM-only 3D velocity dispersion (`Group[].VelDispDM`, mass-weighted `sigma_3D`) exceeds `MinVelDispForFOFSeeding`. |
| `BH_SEED_ON_POTENTIAL_POSITION` | `BLACKHOLE_SEEDING`, `EVALPOTENTIAL` | Changes donor-cell *selection*, not a seeding channel itself: the donor becomes the densest gas cell (searched over all particle types) within `PotentialDonorSearchNSoft` x the DM softening length of the halo's gravitational-potential minimum, instead of the halo's unrestricted densest gas cell. If no gas cell qualifies within range, seeding is deferred to the next FOF pass rather than falling back to the unrestricted donor. |

Channels are OR'd together if more than one is enabled — every enabled channel is evaluated
unconditionally each pass (not short-circuited), so a seed event records every channel that
fired via `FormationChannel` (see below), not just whichever happened to be checked first. A
halo that already hosts a black hole is never reseeded, regardless of channel.

## Runtime parameters (`param.txt`)

Added under `#ifdef HALO_SEEDING` in `src/io/parameters.c`:

- `TimeOfFirstHaloFinding` — simulation time of the first scheduled halo-finding/seeding check.
- `TimeBetweenHaloFinding` — multiplicative factor applied to `All.NextTimeOfHaloFinding` after each check (geometric schedule, matching the pattern used for outputs).
- `MinHaloMassForFOFSeeding` — mass threshold for `BH_SEED_ON_MASS` (only compiled in if that flag is set).
- `ZeroMetallicityThresholdForFOFSeeding` — metallicity threshold for `BH_SEED_ON_ZERO_METALLICITY`.
- `MinVelDispForFOFSeeding` — DM velocity-dispersion threshold for `BH_SEED_ON_VELDISP`.
- `PotentialDonorSearchNSoft` — donor-search cutoff radius for `BH_SEED_ON_POTENTIAL_POSITION`, in units of the DM gravitational softening length. Default suggested range 3-5 (a force-resolution scale, not a halo-boundary one like `LinkL`); **untuned against real halos, worth revisiting**.
- `BlackHoleSeedMass` — mass assigned to newly spawned black hole particles (`BLACKHOLE_SEEDING` only).

## Where it runs in the step loop

Seeding is checked in `src/main/run.c`, inside the main loop, after
`calculate_non_standard_physics_end_of_step()`:

```c
if(All.Time >= All.NextTimeOfHaloFinding && All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
  {
    int num_seed_events = fof_seeding_list(halo_seed_events, MAX_HALO_SEED);
#ifdef BLACKHOLE_SEEDING
    seed_black_holes_from_events(halo_seed_events, num_seed_events);
#endif
    All.NextTimeOfHaloFinding *= All.TimeBetweenHaloFinding;
  }
```

The check requires all particles to be synchronised (`HighestActiveTimeBin == HighestOccupiedTimeBin`)
because the donor gas-cell indices identified during the FOF pass are only valid at that point,
and a new particle can only be safely injected there — before the next domain decomposition.

## Code structure

- **`src/fof/fof.c`** — during FOF group-property accumulation, each group now tracks (under `HALO_SEEDING`):
  - the densest gas cell in the group (candidate seeding donor), including its ID, task, and local index;
  - if `BH_SEED_ON_ZERO_METALLICITY` is set, the most metal-enriched gas cell's metallicity, so a halo is only judged "pristine" once every fragment of the group (across MPI tasks) has been accounted for;
  - if `BH_SEED_ON_VELDISP` is set, a DM-only mass-weighted velocity sum (`Group[].VelDM`) and sum of `|v|^2` (`Group[].VelDispDM`), converted in-place to the DM bulk velocity and 3D velocity dispersion (`sigma_3D = sqrt(mean(v^2) - mean(v)^2)`) in `fof_finish_group_properties()`, following the same accumulate/merge/finalize pipeline as every other group property; `-1` if the group has no DM;
  - if `BH_SEED_ON_POTENTIAL_POSITION` is set, the group's gravitational-potential minimum (`Group[].MinPotential`/`MinPotentialPos`, searched over all particle types, not DM-only) via the same pipeline, converted from `FirstPos`-relative to an absolute periodic-wrapped position at the end.
  - Each gas cell is also tagged with `HostHaloMass`, the FOF mass of its host halo (0 if not in a halo) — refreshed at every FOF pass and usable elsewhere in the code. See `sf_threshold_halo_mass.md` for one such use: raising the star-formation density threshold for gas in low-mass halos (`SF_THRESHOLD_HALO_MASS_DEPENDENT`).

- **`src/fof/fof_seeding.c`** (substantially reworked) — builds the list of halos meeting the enabled seeding channel(s) and produces `HaloSeedEvent` entries (via `fof_seeding_list()`); every enabled channel is evaluated unconditionally (not short-circuited) so `HaloSeedEvent.FormationChannel` records every channel that fired. Calls `fof_seeding_tag_host_halo_mass()` to propagate `HostHaloMass` to gas cells, and — if `BH_SEED_ON_POTENTIAL_POSITION` is set — `fof_seeding_find_potential_donor()`, a new three-phase pass that finds the densest gas cell within the donor-search cutoff of each group's potential minimum:
  1. broadcasts each group's finalized `MinPotentialPos` out to every task holding a fragment of it (mirroring `fof_seeding_tag_host_halo_mass()`'s own query/route/answer/return sequence);
  2. each task searches its own local gas cells within the cutoff and records its own best candidate;
  3. local best candidates are routed back to the group's owning task, which keeps the densest one received.

  Donor selection then reads `Group[].MaxGasDensNearPotential*` instead of `Group[].MaxGasDens*`, deferring (not falling back) if no candidate was found.

- **`src/fof/fof_seeding.h`** (new) — declares `fof_seeding_list()`, `HaloSeedEvent` (including the `FormationChannel` bitmask and the `BH_SEED_CHANNEL_*` constants: `_MASS`, `_ZERO_METALLICITY`, `_VELDISP` — there is no `_POTENTIAL_POSITION` constant, consistent with that flag changing donor *selection* rather than being its own channel), `seed_black_holes_from_events()`, and the `HaloSeedRegistry` type used for restart persistence. `MAX_HALO_SEED` itself is `#define`d in `src/main/run.c`, not here.

- **`src/fof/fof_seeding_registry.c`** (new) — maintains `HaloSeeds`, a registry of already-seeded halos (so they aren't reseeded), and provides `fof_seeding_registry_io()` for restart read/write.

- **`src/blackholes/bh_seed.c`** (new) — `seed_black_holes_from_events()`: consumes the seed-event list and spawns black hole particles (mass-conserving, at the synchronisation point described above), guarded by `#if defined(HALO_SEEDING) && defined(BLACKHOLE_SEEDING)`. `spawn_black_hole_from_cell()` also captures the donor cell's velocity and metallicity (before either is mutated) and populates the `Bh_Particle_Data` formation-history fields below.

- **`src/blackholes/bh.h`** — `Bh_Particle_Data` gains four fields under `#ifdef HALO_SEEDING`: `FormationTime` (`All.Time` at seeding), `FormationMetallicity` (donor cell's metal mass fraction at seeding; `-1` if `METALS` is off), `FormationChannel` (bitmask of `BH_SEED_CHANNEL_*`, copied from the triggering `HaloSeedEvent`), and `DonorVelocity` (donor cell's velocity at seeding). These persist across restarts automatically, since `restart.c` already writes/reads the whole `Bh_Particle_Data` struct as raw bytes.

- **`src/io/restart.c`** — `fof_seeding_registry_io(&HaloSeeds, modus)` persists the seed registry across restarts (alongside, not replacing, the existing black hole particle restart I/O).

- **`src/io/io_fields.c`** — adds `HostHaloMass` (`HHMA`) as an output field on gas cells when `HALO_SEEDING` is enabled, and (under `#if defined(HALO_SEEDING) && defined(BLACKHOLES)`) the four `Bh_Particle_Data` formation-history fields as black-hole-only output fields: `BlackHoleFormationTime` (`BFTI`), `BlackHoleFormationMetallicity` (`BFMZ`), `BlackHoleFormationChannel` (`BFCH`), `BlackHoleDonorVelocity` (`BFDV`).

- **`src/init/begrun.c`** — initialises `All.NextTimeOfHaloFinding` from `All.TimeOfFirstHaloFinding` on fresh starts (not on restart, where it's restored from the restart file).

- **`src/init/init.c`** — `fof_seeding_init(RestartFlag)` is now guarded with `if(RestartFlag != 1)`, mirroring the `All.NextTimeOfHaloFinding` guard in `begrun.c`. Previously this call ran unconditionally and wiped the halo-seed registry that `loadrestart()` had just restored from the restart file, on every `RestartFlag == 1` restart. `RestartFlag == 2` (snapshot continuation) still reseeds the registry from scratch, deliberately — `Group[n].LenType[5] > 0` independently catches any halo whose black hole is currently FoF-linked into it, which any snapshot-restored particle set will show correctly on the next FoF pass.

## Data-integrity fix: BhP[]/SP[] not following P[] during FoF/SubFind output reordering

Every time a snapshot is written with FOF enabled, `fof_subfind_exchange()` in
`src/fof/fof_distribute.c` redistributes `P[]`/`PS[]`/`SphP[]` across MPI tasks to put particles
into group order for output, then reverts the exchange afterward. Until commits `9afdd67` +
`13ae729`, this routine never moved the auxiliary `BhP[]`/`SP[]` arrays (black hole and star
particle data, respectively) alongside a relocated Type-5/Type-4 particle. A relocated BH or
star's `P[].BhID`/`SID` then pointed at the sending task's own local array index — meaningless
on the receiving task — which surfaced as a segfault in `fill_write_buffer()`'s `A_BH` branch
(`src/io/io.c`), since it scans `BhP[]` with no bounds check against `NumBhs`.

This directly affects the trustworthiness of any BH/star snapshot output whenever FOF/SubFind
reordering is active — i.e. essentially always when `HALO_SEEDING`/`BLACKHOLE_SEEDING` are on.
The fix mirrors `domain_exchange()`'s already-working `BhP[]`/`SP[]` handling
(`src/domain/domain_exchange.c`), adapted to `fof_subfind_exchange()`'s per-type do-while
structure: `P[].BhID`/`SID` are resynced at the top of each type's cycle (since earlier types'
bulk `memmove()`s in the same function can silently relocate a BH/star locally without anything
else knowing), the export/import paths copy the auxiliary entries alongside the main particle
buffer with swap-with-last reclamation of vacated slots, `BhP[]`/`SP[]` are grown via
`reallocate_memory_maxpartbhs()`/`reallocate_memory_maxpartstars()` as needed, and a final
comprehensive `PID` resync pass runs at the end of the function. A companion fix in
`src/utils/allocate.c` makes those two reallocation routines `memset()` newly-grown memory,
matching `allocate_memory()`'s existing behaviour on initial allocation. Verified against a
synthetic reproduction (segfault reproduced 3/3 trials pre-fix; fix confirmed across 2/3/4 MPI
ranks and multiple snapshots with all seeded BH IDs intact).

## Known limitations / not yet validated

- **`fof_seeding_find_potential_donor()`'s gather-back step** (phase 3 above) is the one piece
  of this feature that's genuinely new MPI communication rather than a close mirror of
  something already working elsewhere in the file. Verified only by compiling successfully
  (including with all four seeding channels enabled together, and confirming the
  `EVALPOTENTIAL` `#error` guard fires when it's missing) — not run against a real halo
  catalogue. Recommend testing it in isolation (e.g. a small box with a known,
  hand-checkable gas configuration near a potential minimum) before trusting it in a full run.
  This function was also the site of a LIFO-allocator-ordering crash (`0491776`): `myfree(query)`
  was originally called before the loop consuming `candidates` had finished, but `query` was
  allocated *after* `candidates` on Arepo's strict LIFO `mymalloc`/`myfree` stack, so it had to be
  freed first — fixed by reordering the two calls.
- **`BH_SEED_ON_ZERO_METALLICITY` was silently broken until `c0f1919`**: `bh_seed.c` and `fof.c`'s
  zero-metallicity block referenced a non-existent `SphP[].Metals` field instead of the
  `GasMetals` macro (`PConservedScalars[METALS_INDEX]`). This went uncaught because the channel
  is off in `CosmoConfig.sh`; it has since been fixed and build-verified with the channel enabled,
  but has no run history yet — treat it the same as the velocity-dispersion/potential-position
  channels for validation purposes.
- **`PotentialDonorSearchNSoft` default (5.0)** is an untuned starting point from the design
  discussion's suggested 3-5 range, not derived from testing against real halos.
- **Registry/`LenType[5]` joint failure**: a black hole that both (a) has drifted to a
  no-longer-registered `MinID` due to a marginal merge/split, and (b) isn't currently
  FoF-linked into its own host group — narrow, unverified frequency, accepted rather than
  hardened against.
- **`R_vir`-based donor-search cutoff**: the softening-length-based cutoff
  (`PotentialDonorSearchNSoft` x DM softening) is a deliberate placeholder for a more
  physically principled, halo-scaled version. Deferred because `R_vir` isn't computed in this
  lightweight, SubFind-independent seeding pass; revisit only if the softening-based cutoff
  proves too coarse in practice.
- **Ongoing post-formation repositioning** toward the potential minimum (not just at the moment
  of seeding) is established elsewhere in the literature but deliberately not implemented here;
  `DonorVelocity` is preserved on each seeded black hole specifically to allow checking, after
  the fact, whether donor-vs-bulk velocity at formation actually mattered.
- **General descriptor-based seeding framework** (separating "what is the physical state of
  this candidate" from "should a black hole form here" as a formal, pluggable abstraction):
  explicitly out of scope until there are 3+ *validated* seeding criteria, not merely
  implemented ones — mass and zero-metallicity seeding have some run history, but
  velocity-dispersion and potential-anchored donor selection do not yet.

## Post-merge robustness fixes (production crash investigation, August 2026)

Two crashes surfaced in real Setonix production runs after the merge, both in the registry
growth path (`fof_seeding_registry.c`/`fof_seeding.c`) rather than in feature logic. Neither
was reachable in normal short test runs — both needed a restart with zero pre-existing
allocation headroom (see below) to trigger.

- **`halo_seed_registry_grow()` reallocated the wrong pointer.** It called
  `myrealloc_movable(&r->ids, ...)` — passing the *address of the struct field* — instead of
  `myrealloc_movable(r->ids, ...)`, which takes the current pointer *value*.
  `myrealloc_movable()`/`myfree_movable()` look up a block by the pointer value tracked
  internally, so passing the wrong thing corrupted the memory-arena bookkeeping the next time
  anything tried to free or resize that block. Fixed by passing `r->ids` and reassigning the
  result, matching every other `myrealloc_movable()` call site in the codebase.
- **Non-movable scratch alive behind `HaloSeeds.ids` during registry growth.** `mymalloc`/
  `myfree` is a strict LIFO stack: a `movable` block (like `HaloSeeds.ids`) can only grow via
  `myrealloc_movable()` if every block allocated *after* it on the arena is also movable — any
  `mymalloc()` (non-movable) block still alive above it blocks the resize. `fof_seeding_list()`
  was calling `halo_mark_seeded()` (which can trigger registry growth) *before* freeing its own
  non-movable scratch buffers (`bdispls`/`bcounts`/`counts`/`local_events`) and non-movable-ish
  `FOF_GList`/`FOF_PList`/`Group`/`PS`. Fixed by reordering so all of that cleanup happens
  before `halo_mark_seeded()` runs. This only manifested once enough halos had been seeded
  across a run's history to exhaust the registry's initial headroom and force a real growth —
  hence only appearing well into long production runs, never in short test runs.

Both fixes are on `merge-fof-into-star-feedback` only for now (commits `43ed2f5`, `d1677cb`) --
`fof_seeding_registry.c` doesn't exist on `Star_feedback_radiation` yet, so there was nothing to
cherry-pick there; port these alongside the wider `HALO_SEEDING` feature merge, not before.

A separate, related fix landed in `domain_rearrange_particle_sequence()`
(`src/domain/domain_rearrange.c`, commit `80fe435`): when a particle is eliminated under
`REFINEMENT_MERGE_CELLS` and it is itself a black hole (e.g. a BH merged away by `bh_merger()`),
its own `BhP[]` slot was never reclaimed, leaving `NumBhs` permanently out of sync with the
true live-BH count. This affects any black hole regardless of how it was created, seeded ones
included. This fix *is* portable and has its own PR against `Star_feedback_radiation`
(see the general fix-porting workflow — not seeding-specific, so not detailed further here).

## Further robustness fixes (mesh-connectivity + fossilized BhP[] crashes, August 2026)

A second wave of production crashes surfaced after the fixes above, none of them specific to
this feature's own logic — they trace to two separate, unrelated root causes that this feature
happens to expose more often than most (frequent on-the-fly FOF passes each force their own
`domain_Decomposition()`, and BH-driven derefinement of donor cells adds extra mesh churn).

**Root cause 1: fossilized zombie `BhP[]` entries from pre-event-based seeding.** Git archaeology
(`git log --follow -- src/blackholes/bh_seed.c`) found that before commit `7d15250` ("Event-based
BH seeding"), an earlier "basic form" of seeding (commit `4347420`) had its BH-structure
initialisation block entirely commented out — a gas cell was flipped to `Type=5` but `NumBhs` was
never incremented and `P[].BhID` was never assigned, defaulting to `0`. Every BH seeded under that
old code collided on `BhID=0`/`BhP[0]`, leaving old restart checkpoints with a permanently
under-initialised slot (`TimeBinBh`, `SofteningType`, `Hsml` all garbage) that later code paths
assumed was always valid. This is fossilized checkpoint history, not a live bug — a dedicated
BH-seeding stress test against current `bh_seed.c` (many seed/regrow events, 1- and 2-task runs)
found nothing wrong with the seeding code itself. Fixed by self-healing wherever this garbage
surfaces, loudly, instead of crashing:

- `src/blackholes/bh_update.c` (`bh_reconstruct_timebins()`, commit `c135099`): `TimeBinBh` outside
  `[0, TIMEBINS)` is reset to bin `0` (not skipped — an earlier skip-based attempt permanently
  orphaned the BH from timebin scheduling forever).
- `src/blackholes/bh_density.c` (commit `7908521`, hardened `98a36d8`): degenerate
  `SofteningType` is reset to the BH type default. The first version trusted `PPB(i).Type` for the
  fallback lookup, which is itself garbage for these entries (`Type` is `unsigned char`, so values
  up to 255 are possible and index far past `SofteningTypeOfPartType[NTYPES]`); hardened to use a
  hardcoded BH-type literal (`5`, matching `bh_seed.c`'s convention) instead.
- `src/blackholes/bh_accretion.c` (`bh_accretion_rate()`, commit `fc29fef`): a Bondi-rate
  denominator that's `<=0` (or NaN) only happens when the local gas has zero sound speed and zero
  velocity, or `GasInternalEnergy` is corrupted negative — the same zombie-entry signature. Treated
  as "no accretion" (matching the existing `GasDensity<=0` branch immediately above) instead of
  terminating.

**Root cause 2: a rare, latent bug in AREPO's own mesh-connectivity exchange.** Confirmed via
`git log --follow` that `src/domain/domain_DC_update.c` and
`src/mesh/voronoi/voronoi_dynamic_update.c` have never been modified since the original public
AREPO import — this is genuine upstream behaviour, not a Solas regression. `domain_exchange_and_update_DC()`
rebuilds the entire `DC[]` connection array on every domain decomposition; under certain edge
cases a connection can fail to cleanly reattach to any cell's chain during that rebuild, without
being returned to the free list either (the free list is rebuilt from the *tail* of the array
before the per-particle reattachment loop runs) — it becomes an inert "dead zone" slot that no
cell references, and can sit dormant for tens of thousands of sync-points before some unrelated
cell's connection-chain walk coincidentally routes through that same index, at which point it
looks like a self-loop or a foreign-particle collision. Three concrete edge cases were found and
fixed, all in `src/domain/domain_DC_update.c`:

- `domain_exchange_and_update_DC()`'s first exchange round (transcribing task/index) was missing
  the `DC[i].index >= 0` guard the second round already had, so an already-removed connection
  could be re-sent and crash the receiver's `trans_table[]` lookup (commit `481f3d6`).
- `domain_mark_in_trans_table()`'s per-cell chain walk can hit a slot whose `.next` points back to
  itself (`q == DC[q].next`) — confirmed via a targeted diagnostic to be a genuinely inactive,
  otherwise-normal cell whose claimed connection range has been silently reallocated to an
  unrelated particle partway through. Root mechanism not fully pinned down; mitigated by breaking
  out of the chain walk at that point with a loud warning instead of terminating (commit `9b6fa3f`).
- That break leaves `DC[q].next` unset for the rest of the chain, which two downstream loops read
  as a destination-task scratch value without validating it — first
  `domain_exchange_and_update_DC()`'s second exchange round (a stale value `>= NTask` crashed
  `terminate()`, commit `4c5fe31`), then its final ID-based merge-join reattaching connections to
  local particles (a stale/orphaned entry with no matching local ID hit `terminate("strange")`,
  commit `186deaa`). Both now drop the offending connection with a warning instead of crashing.

**Investigated and ruled out** as the missing-cleanup culprit: `spawn_black_hole_from_cell()`
(`bh_seed.c`) never removes a gas cell — the donor persists with reduced mass — so it was never a
candidate for a missing `voronoi_remove_connection()` call. Both of the two paths that *do* fully
remove a gas cell already call it correctly: `voronoi_derefinement.c:711` and
`starformation.c`'s `convert_cell_into_star()` (line ~346). A BH-free, pure-hydro synthetic
reproduction (no gravity, no `HALO_SEEDING`/`BLACKHOLES`/`FOF`, just `REFINEMENT_SPLIT_CELLS`/
`MERGE_CELLS` driven by a standing compression wave, genuine sustained refine/derefine churn over
a full run) came through clean, suggesting the trigger needs either self-gravity's cross-task
particle migration pattern or far more exposure than that test achieved — not yet confirmed either
way.

**What increases exposure, and what doesn't fix it:** the bug's frequency scales with how often
`domain_Decomposition()` runs — `ActivePartFracForNewDomainDecomp` (currently `0.01` in
`param.txt`, an aggressive/frequent setting) and, independently, every on-the-fly FOF/seeding pass
(`src/fof/fof.c` calls `domain_Decomposition()` unconditionally, outside that threshold). Relaxing
`ActivePartFracForNewDomainDecomp` to something like `0.1` would reduce exposure at no physics
cost; reducing `TimeBetweenHaloFinding`'s frequency would also reduce exposure but has a real
physics cost (delayed BH seeding) and shouldn't be tuned purely to dodge this. Two production runs
with identical IC/`param.txt`/compiler flags and even matched task count (Setonix, 8 tasks, crashed
repeatedly vs. OzSTAR, 8 and 16 tasks, both clean) show the trigger is sensitive to something below
that level of reproducibility — most plausibly floating-point non-associativity across different
machines/compilers/task counts compounding, via this being a chaotically-sensitive self-gravitating
system, into different halo-assembly timing and therefore different exposure. Not confirmed by a
direct bit-level comparison.

**Status as of August 2026**: all of the above fixes are on `merge-fof-into-star-feedback`. Two are
also on branches with **open, unmerged** PRs against `Star_feedback_radiation` — PR #34
(`d8b2742`/superseded-by-`c135099`, `7908521` — **not yet updated with the `98a36d8` hardening**)
and PR #35 (`481f3d6` only). The three later `domain_DC_update.c` self-heals (`9b6fa3f`,
`4c5fe31`, `186deaa`) and `bh_accretion.c`'s fix (`fc29fef`) have not been ported to
`Star_feedback_radiation` at all yet — production runs on that branch would still hit the original
crashes. The `ActivePartFracForNewDomainDecomp` relaxation has not been committed to the tracked
`param.txt` either.

## Development history (chronological, oldest first)

The feature evolved through several rounds of bug-fixing before settling into its current form:

1. Initial FOF-seeding pass produced a basic version that changed a gas cell's type in place rather than spawning a real black hole particle.
2. Seeding logic was split out of `bh_seed.c`/`fof_seeding.c` into a cleaner interface, with the new registry (`fof_seeding_registry.c`) added for tracking already-seeded halos.
3. A series of fixes addressed: an uninitialised `NextTimeOfHaloFinding`, mesh reconstruction/MPI crashes during seeding, ordering bugs when freeing `Group`/`PS` arrays, a memory leak in registry restart reads, and a crash caused by the registry I/O routine freeing the live registry during restart writes.
4. Spawning was reworked to be event-based and mass-conserving, and moved to the natural synchronisation point in the step loop (mirroring star particle spawning).
5. `HostHaloMass` tagging was added so gas cells carry their host halo's FOF mass, fixing an initial bug where it was always zero (secondary particle types weren't being linked correctly in the on-the-fly FOF pass).
6. Selectable seeding channels (`BH_SEED_ON_MASS`, `BH_SEED_ON_ZERO_METALLICITY`) were added as the final piece, replacing a single hard-coded criterion.
7. `FIND_HALOS` was renamed to `HALO_SEEDING` partway through this history, and `BLACKHOLE_SEEDING` was introduced as a separate flag from the seeding-detection logic itself.

## Note on the merge into `Star_feedback_radiation`

`Star_feedback_radiation` had independently developed a similar but simpler mechanism under
`FIND_HALOS`: a periodically scheduled call to `fof_seeding()` from
`calculate_non_standard_physics_prior_mesh_construction()`, using the same
`TimeOfFirstHaloFinding` / `NextTimeOfHaloFinding` / `TimeBetweenHaloFinding` field names but
without the event registry, seeding channels, or proper black hole spawning.

During the merge, `FIND_HALOS` was retired in favour of this branch's `HALO_SEEDING` /
`BLACKHOLE_SEEDING` implementation:

- The old scheduled call in `run.c` was removed; the `fof`-branch call at the post-end-of-step
  synchronisation point (see above) is the one now used.
- Duplicate struct members and object-file rules that `Star_feedback_radiation` had already
  evolved independently (e.g. its own more advanced black hole physics call sequence, STARS
  object list, and library-linking configuration) were kept as-is rather than overwritten,
  since the `fof` branch's versions of those particular pieces predated `Star_feedback_radiation`'s
  own development and were superseded by it.
- `CosmoConfig` (deleted on `Star_feedback_radiation` in favour of `CosmoConfig.sh`) was kept
  deleted; the `fof` branch's edits to the old file were superseded.

No functional code from either branch's genuinely distinct physics (star feedback/radiation on
one side, halo/BH seeding on the other) was dropped — only duplicated or superseded pieces of
build configuration and scheduling logic were consolidated.
