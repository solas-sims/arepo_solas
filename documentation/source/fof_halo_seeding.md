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
  - Each gas cell is also tagged with `HostHaloMass`, the FOF mass of its host halo (0 if not in a halo) — refreshed at every FOF pass and usable elsewhere in the code, e.g. to select star formation mode by halo mass.

- **`src/fof/fof_seeding.c`** (substantially reworked) — builds the list of halos meeting the enabled seeding channel(s) and produces `HaloSeedEvent` entries (via `fof_seeding_list()`); every enabled channel is evaluated unconditionally (not short-circuited) so `HaloSeedEvent.FormationChannel` records every channel that fired. Calls `fof_seeding_tag_host_halo_mass()` to propagate `HostHaloMass` to gas cells, and — if `BH_SEED_ON_POTENTIAL_POSITION` is set — `fof_seeding_find_potential_donor()`, a new three-phase pass that finds the densest gas cell within the donor-search cutoff of each group's potential minimum:
  1. broadcasts each group's finalized `MinPotentialPos` out to every task holding a fragment of it (mirroring `fof_seeding_tag_host_halo_mass()`'s own query/route/answer/return sequence);
  2. each task searches its own local gas cells within the cutoff and records its own best candidate;
  3. local best candidates are routed back to the group's owning task, which keeps the densest one received.

  Donor selection then reads `Group[].MaxGasDensNearPotential*` instead of `Group[].MaxGasDens*`, deferring (not falling back) if no candidate was found.

- **`src/fof/fof_seeding.h`** (new) — declares `fof_seeding_list()`, `HaloSeedEvent` (including the `FormationChannel` bitmask and the `BH_SEED_CHANNEL_*` constants), `MAX_HALO_SEED`, and the `HaloSeedRegistry` type used for restart persistence.

- **`src/fof/fof_seeding_registry.c`** (new) — maintains `HaloSeeds`, a registry of already-seeded halos (so they aren't reseeded), and provides `fof_seeding_registry_io()` for restart read/write.

- **`src/blackholes/bh_seed.c`** (new) — `seed_black_holes_from_events()`: consumes the seed-event list and spawns black hole particles (mass-conserving, at the synchronisation point described above), guarded by `#if defined(HALO_SEEDING) && defined(BLACKHOLE_SEEDING)`. `spawn_black_hole_from_cell()` also captures the donor cell's velocity and metallicity (before either is mutated) and populates the `Bh_Particle_Data` formation-history fields below.

- **`src/blackholes/bh.h`** — `Bh_Particle_Data` gains four fields under `#ifdef HALO_SEEDING`: `FormationTime` (`All.Time` at seeding), `FormationMetallicity` (donor cell's metal mass fraction at seeding; `-1` if `METALS` is off), `FormationChannel` (bitmask of `BH_SEED_CHANNEL_*`, copied from the triggering `HaloSeedEvent`), and `DonorVelocity` (donor cell's velocity at seeding). These persist across restarts automatically, since `restart.c` already writes/reads the whole `Bh_Particle_Data` struct as raw bytes.

- **`src/io/restart.c`** — `fof_seeding_registry_io(&HaloSeeds, modus)` persists the seed registry across restarts (alongside, not replacing, the existing black hole particle restart I/O).

- **`src/io/io_fields.c`** — adds `HostHaloMass` (`HHMA`) as an output field on gas cells when `HALO_SEEDING` is enabled, and (under `#if defined(HALO_SEEDING) && defined(BLACKHOLES)`) the four `Bh_Particle_Data` formation-history fields as black-hole-only output fields: `BlackHoleFormationTime` (`BFTI`), `BlackHoleFormationMetallicity` (`BFMZ`), `BlackHoleFormationChannel` (`BFCH`), `BlackHoleDonorVelocity` (`BFDV`).

- **`src/init/begrun.c`** — initialises `All.NextTimeOfHaloFinding` from `All.TimeOfFirstHaloFinding` on fresh starts (not on restart, where it's restored from the restart file).

- **`src/init/init.c`** — `fof_seeding_init(RestartFlag)` is now guarded with `if(RestartFlag != 1)`, mirroring the `All.NextTimeOfHaloFinding` guard in `begrun.c`. Previously this call ran unconditionally and wiped the halo-seed registry that `loadrestart()` had just restored from the restart file, on every `RestartFlag == 1` restart. `RestartFlag == 2` (snapshot continuation) still reseeds the registry from scratch, deliberately — `Group[n].LenType[5] > 0` independently catches any halo whose black hole is currently FoF-linked into it, which any snapshot-restored particle set will show correctly on the next FoF pass.

## Known limitations / not yet validated

- **`fof_seeding_find_potential_donor()`'s gather-back step** (phase 3 above) is the one piece
  of this feature that's genuinely new MPI communication rather than a close mirror of
  something already working elsewhere in the file. Verified only by compiling successfully
  (including with all four seeding channels enabled together, and confirming the
  `EVALPOTENTIAL` `#error` guard fires when it's missing) — not run against a real halo
  catalogue. Recommend testing it in isolation (e.g. a small box with a known,
  hand-checkable gas configuration near a potential minimum) before trusting it in a full run.
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
