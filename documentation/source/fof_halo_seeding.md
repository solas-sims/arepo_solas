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
| `HALO_SEEDING` | `FOF` | Enables tracking of halo properties needed for seeding (densest gas cell, host halo mass, optionally metallicity) during each FOF pass, and the seed-event list mechanism. |
| `BLACKHOLE_SEEDING` | `HALO_SEEDING` | Enables spawning black hole particles from seed events. Requires at least one seeding channel below. |
| `BH_SEED_ON_MASS` | `BLACKHOLE_SEEDING` | Seeding channel: seed once a halo's FOF mass exceeds `MinHaloMassForFOFSeeding`. |
| `BH_SEED_ON_ZERO_METALLICITY` | `BLACKHOLE_SEEDING`, `METALS` | Seeding channel: seed "pristine" halos, i.e. halos whose most metal-enriched gas cell still has a metal mass fraction at or below `ZeroMetallicityThresholdForFOFSeeding`. |

Channels are OR'd together if more than one is enabled. A halo that already hosts a black hole
is never reseeded, regardless of channel.

## Runtime parameters (`param.txt`)

Added under `#ifdef HALO_SEEDING` in `src/io/parameters.c`:

- `TimeOfFirstHaloFinding` — simulation time of the first scheduled halo-finding/seeding check.
- `TimeBetweenHaloFinding` — multiplicative factor applied to `All.NextTimeOfHaloFinding` after each check (geometric schedule, matching the pattern used for outputs).
- `MinHaloMassForFOFSeeding` — mass threshold for `BH_SEED_ON_MASS` (only compiled in if that flag is set).
- `ZeroMetallicityThresholdForFOFSeeding` — metallicity threshold for `BH_SEED_ON_ZERO_METALLICITY`.
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
  - if `BH_SEED_ON_ZERO_METALLICITY` is set, the most metal-enriched gas cell's metallicity, so a halo is only judged "pristine" once every fragment of the group (across MPI tasks) has been accounted for.
  - Each gas cell is also tagged with `HostHaloMass`, the FOF mass of its host halo (0 if not in a halo) — refreshed at every FOF pass and usable elsewhere in the code, e.g. to select star formation mode by halo mass.

- **`src/fof/fof_seeding.c`** (substantially reworked, 405 lines changed) — builds the list of halos meeting the enabled seeding channel(s) and produces `HaloSeedEvent` entries (via `fof_seeding_list()`), and calls `fof_seeding_tag_host_halo_mass()` to propagate `HostHaloMass` to gas cells.

- **`src/fof/fof_seeding.h`** (new) — declares `fof_seeding_list()`, `HaloSeedEvent`, `MAX_HALO_SEED`, and the `HaloSeedRegistry` type used for restart persistence.

- **`src/fof/fof_seeding_registry.c`** (new) — maintains `HaloSeeds`, a registry of already-seeded halos (so they aren't reseeded), and provides `fof_seeding_registry_io()` for restart read/write.

- **`src/blackholes/bh_seed.c`** (new, 219 lines) — `seed_black_holes_from_events()`: consumes the seed-event list and spawns black hole particles (mass-conserving, at the synchronisation point described above), guarded by `#if defined(HALO_SEEDING) && defined(BLACKHOLE_SEEDING)`.

- **`src/io/restart.c`** — `fof_seeding_registry_io(&HaloSeeds, modus)` persists the seed registry across restarts (alongside, not replacing, the existing black hole particle restart I/O).

- **`src/io/io_fields.c`** — adds `HostHaloMass` (`HHMA`) as an output field on gas cells when `HALO_SEEDING` is enabled.

- **`src/init/begrun.c`** — initialises `All.NextTimeOfHaloFinding` from `All.TimeOfFirstHaloFinding` on fresh starts (not on restart, where it's restored from the restart file).

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
