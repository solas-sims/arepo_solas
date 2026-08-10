# Technical Reference: `src/stars/` and `src/blackholes/`

**Branch documented:** `merge-fof-into-star-feedback`

**Last updated:** 6 August 2026, against commit `fb06d35`. This supersedes the project wiki's
[Technical-Reference](https://github.com/solas-sims/arepo_solas/wiki/Technical-Reference) page,
which is pinned to commit `f692626` (24 July 2026) — this page picks up the ~67 commits since
then. The wiki page now points here for the current version rather than duplicating it; see the
note at the bottom of this page.

**Scope:** File-by-file, function-level detail on the two new physics modules (`src/stars/`,
`src/blackholes/`), plus the supporting star-formation and infrastructure changes referenced
there. For the on-the-fly FOF halo-seeding mechanism itself (event list, seeding channels,
registry), see [fof_halo_seeding.md](fof_halo_seeding.md) — this page cross-references it rather
than repeating it, except where the black-hole side of seeding (`bh_seed.c`) is concerned.

This is the "what does this file actually do" document. Read it with the code open — line
numbers are given so you can jump straight to a function. All line numbers below are current
(HEAD, `fb06d35`), not diff line numbers.

---

## 1. `src/stars/` — star particle framework

### 1.1 Data model (`star.h`, `star_particle.h`, `star_tables.h`)

Every simulation particle of the star type gets a companion struct `Star_Particle_Data` (in
`star.h`), accessed via the `SPP(i)` / `PPS(i)` macros. Key fields, each gated by a `#ifdef` so
they only exist in builds that need them:

- `MassOfStar`, `PhysicalAge_yr`, `Hsml`, `NgbsMass`, `NgbsVolume` — basic bookkeeping and the
  SPH-like smoothing length/neighbor mass used to find a star's host gas.
- `Active` / `WithFeedback` — lifecycle flags: `Active` is a permanent state (`STAR_UNBORN` = 0,
  `STAR_ACTIVE` = 1, `STAR_INACTIVE` = −1), while `WithFeedback` marks per-timestep whether the
  star produced any feedback this step. Feedback stages are tracked separately as `STAR_MS` /
  `STAR_SN` / `STAR_POST_SN` in `Star_Feedback.Stage`.
- `Birthtime` (line 217, added in `52d28b2`) — `All.Time` at the step the star activated
  (`Active` transitions from `STAR_UNBORN`). See §1.6 for how this replaced incremental age
  accumulation.
- `MechanicalFeedback` (a nested `Mechanical_Feedback` struct) — the per-star payload of mass,
  metals, H/He, momentum, and radiated energy/photons **by waveband** delivered to gas each step.
- `MassToDrain` — only present under `INDIVIDUAL_STAR_BY_STAR_FORMATION`; tracks how much gas
  mass still needs to be siphoned out of surrounding cells to "pay for" a star that already
  formed (see §3).
- `TimeToSN` / `NextSNEnergy` — only present with `TREE_BASED_TIMESTEPS` + `SUPERNOVAE`; renamed
  from `TimeSN_yr` in `03ab10f` (see §1.4) — the precomputed *time until* this star's SN (not the
  absolute lifetime) and its already-unit-converted injection energy, so the tree-based timestep
  criterion can see the SN coming and shrink the timestep in advance.

Two mass thresholds are defined in `star.h` (both in M☉): `LOWEST_MASS_FEEDBACK = 2` (below this
a star contributes no feedback at all) and `LOWEST_MASS_SN = 8` (below this a star never explodes
as a supernova). `SN_ENERGY` (`star.h:36`, added in `18a87f2`) is now the single `#define
1.0e51` used throughout `star_interpolation.c` and `star_feedback.c`, replacing scattered literal
`1e51`/`1.0e51` occurrences.

Three companion structs handle **delivering** feedback:
- `Mechanical_Feedback` — the physical payload: position/velocity, wind mass/metals/momentum
  (`WINDS`), per-waveband radiated energy and photons as `WavebandData` pairs
  (`STAR_RADIATION_ACTIVE`), and SN ejecta mass/metals/energy (`SUPERNOVAE`). **Struct shape
  changed** in the H/He injection work (§1.7 below): `HLoss`/`HeLoss` (wind) and
  `SN_HLoss`/`SN_HeLoss` (SN) fields were inserted between the existing mass and metals fields.
  Because this struct is broadcast wholesale in `star_density.c`'s `data_in`/`data_out`
  (`in->MechanicalFeedback = SP[i].MechanicalFeedback;`), the change altered the wire size of
  every `Mechanical_Feedback_Data` MPI message — all code paths were updated together, but any
  external tooling that deserializes this struct's binary layout needs to match.
- `Mechanical_Feedback_Data` — wraps a payload with the star's and host cell's local/MPI-task
  indices, so a feedback event can be routed to whichever process owns the target gas cell.
- `Mechanical_Feedback_Events` — a dynamically-resized event list, managed by
  `feedback_init`/`feedback_allocate`/`feedback_reallocate`/`feedback_free` (`star_update.c:63–85`).
  The global instance `MechanicalFeedbackEvents` lives in `star.c`.

`star_tables.h`/`star_tables.c` define 2D/3D lookup tables (`Age`, `Radius`, `Temperature`,
`MassLossRate`, `SN_MassLoss`, `Flux[WAVEBANDS]`, etc.), indexed by metallicity bin and mass bin,
loaded from an HDF5 file (`load_star_tables()`, `star_tables.c:210`) named by the
`StarTablesFile` parameter — this fork's stand-in for a stellar population synthesis lookup.

### 1.2 Two star-particle "modes" (`STAR_PARTICLES`)

Unchanged since the wiki's last snapshot: `STAR_PARTICLES` = 0/1 selects population-sampled star
particles (IMF-binned across `NBINS = 114` mass bins); `STAR_PARTICLES == 2` selects individual,
resolved stars (§3).

### 1.3 Density / host-finding (`star_density.c`)

Same two-phase tree/ngb-loop pattern as before (`star_density_evaluate1`/`_evaluate2`, driven by
`star_density()`). Two changes since the wiki snapshot:

- `feedback_compare()` (lines 32–44, from `52d28b2`) dropped `HostTask` as the primary sort key —
  now sorts purely by `HostIndex` then `StarTask`/`StarIndex`.
- The Hsml-search and pass-2 re-activation loops (lines ~292–299, ~358–367, from `8a43fe7`) now
  `continue` (skip) over stars with `WithFeedback == 0`, so deactivated/no-feedback stars are no
  longer redundantly walked by the neighbour-search machinery. A new post-pass-2 cleanup loop
  (lines 371–378) `memset()`s every active star's `MechanicalFeedback` payload before the next
  step's `star_prep()` repopulates it, clearing stale data between steps.

### 1.4 Feedback delivery (`star_feedback.c`) — mesh path with a host-only fallback

Delivery-mode selection (`MESH`/`HOST` per event) and the core mesh-face and host-fallback
mechanics are unchanged from the wiki's description — see that page's §1.4 for the full
Kim & Ostriker (2015) cooling-radius / Kimm et al. (2015) terminal-momentum treatment, still
current. What's changed since `f692626`:

- **H/He injection** (§1.7 below) is now threaded through the same kick machinery:
  `struct Feedback_Kick` (lines 18–46) gained `DeltaChem[GRACKLE_SPECIES_NUMBER]` and
  `SN_DeltaChem[GRACKLE_SPECIES_NUMBER]` (generalized chemistry-species delta arrays, replacing
  what would otherwise be separate H/He scalars), indexed by `CHEM_HI/HII/HeI/HeII/HeIII/...`
  constants (`src/main/allvars.h:540–556`). `apply_kick()` (lines 49–86) now loops over all
  Grackle species depositing into a new per-cell buffer `SphP[j].StarChemFeed[s]`
  (`allvars.h:1791`), mirroring how metals are handled.
- **`SN_ENERGY` macro adoption**: lines 143, 199 now compute `E51 = E_SN *
  All.cf_UnitEnergy_in_cgs / SN_ENERGY` instead of dividing by a bare `1e51` literal.
- Wind H/He is injected **neutral** (`CHEM_HI`/`CHEM_HeI`, lines ~718–721); SN H/He is injected
  **ionized**, on top of the swept-host-cell chemistry contribution (`CHEM_HII`/`CHEM_HeIII`,
  lines ~744–748) — a physically deliberate distinction between the two channels.

### 1.5 Stellar evolution interpolation (`star_interpolation.c`) — now log-log

This file changed the most of any single file in the whole diff (302 lines) because two
unrelated things landed here together: the H/He injection plumbing (§1.7) and a restructuring of
*how* the interpolation itself is done.

**Log-space interpolation (`b804bc5`, `4cdede2`, `ecb33d4`).** Previously, mass/metallicity
bracketing used the linear `M_VALUES[]`/`Z_VALUES[]` tables directly. Now two new global tables,
`logZ_VALUES`/`logM_VALUES` (declared `star_tables.h:13–14`), are precomputed once at load time
as `log10()` of the linear tables (`star_tables.c:247`, `:256`, with monotonicity/positivity
guards that `terminate()` at startup if violated) and broadcast to all ranks via `MPI_Bcast`
(`star_tables.c:278–279` — this broadcast is what `4cdede2` added; without it, non-root ranks
would interpolate against zero-initialized log tables). Every bracket-and-interpolate function
now works in log-x space:

- `star_lifetime()`/`lifetime()` (lines 39–75, 77–102) — bracket on `log10(mass)`/`log10(Z)`
  against `logM_VALUES`/`logZ_VALUES`, then interpolate `log10(lifetime)` linearly and transform
  back via `pow(10, ...)` — true log-log (power-law) interpolation, appropriate since stellar
  lifetimes span many orders of magnitude across the mass range.
- `interpolate_mass()`/`interpolate_metallicity()` (lines 233–281, 284–332) and
  `SN_interpolate_mass()`/`SN_interpolate_metallicity()` (lines 337–435, 438–486) — bracket on
  log-x but interpolate the physical quantities themselves (Radius, Temperature, mass-loss rates,
  H/He/MetalsLossRate, WindVelocity, Flux, SN yields) linearly in log-mass/log-metallicity
  *position*, not log-transforming the y-values.

**FSN/DBH clamping fix (`2e921de`).** "FSN"/"DBH" = Failed SN / Direct-Collapse Black Hole —
stellar masses where `SN_MassLoss` is tabulated as 0 because the star collapses without a
genuine explosion. Before this fix, interpolating between a real-SN grid point and a failed-SN/
DCBH grid point silently blended a real yield with zero, producing nonsense. Now, when only one
bracketing endpoint is a real SN, the code **clamps to the nearer grid point in log-distance**
instead of interpolating (lines 401–428 for mass, 471–478 for metallicity).

**`next_SN_time()`** (lines 104–116, `TREE_BASED_TIMESTEPS && SUPERNOVAE`) simplified to always
call `SN_interpolate_metallicity()` and return `tau - a` (time *until* SN) if `SN_MassLoss > 0`,
else `MAX_REAL_NUMBER` for a failed SN/DCBH — the redundant mass-threshold guard that used to be
duplicated both inside and outside this function was consolidated to the single call site in
`star_feedback_compute()` (line 570).

**SN-vs-MS branch reordering (`52d28b2`).** `star_feedback_compute()` (lines 489–581) now checks
the `STAR_SN` stage (`tau < a + dt`) *before* the `STAR_MS`/`else` branch (previously the SN
branch was only reachable via `else`, which didn't correctly account for stars whose lifetime
falls within the *current* timestep). The post-SN check itself changed from `a >= tau + dt` to
`a > tau` (line 505), an off-by-one fix.

**H/He propagation** — `interpolate_age/mass/metallicity()`, `SN_interpolate_mass/metallicity()`,
`star_feedback_compute()`, `units_for_feedback()`, and `star_particle_feedback()` all now
propagate `HLoss`/`HeLoss` (wind) and `SN_HLoss`/`SN_HeLoss` alongside the pre-existing mass/
metals quantities at every stage — see §1.7 for the full data-flow trace.

### 1.6 Timesteps and bookkeeping (`star_update.c`)

- `star_timestep()`, `star_update_timesteps()`, `star_reconstruct_timebins()`,
  `star_update_list_of_active_particles()` — unchanged in role since the wiki snapshot.
- **`deactivate_star(int i)`** (new, `static inline`, lines 205–221, from `8a43fe7`) — consolidates
  what used to be two scattered lines (`Active = STAR_INACTIVE; WithFeedback = 0;`) into a single
  routine that also resets `TimeToSN = MAX_REAL_NUMBER`, `NextSNEnergy = 0.0`, `Hsml = 0.0`,
  `DensityFlag = -1`, `NgbsMass = 0.0`, `NgbsVolume = 0.0`, `HostHydroBin = TIMEBINS`,
  `TimeBinStar = TIMEBINS` — fully parking the star out of the active neighbour-search/timestep
  machinery rather than just flipping a flag. Called from `star_prep()` (line 271) when
  `StarFeedback.State == -1`.
- **`star_prep()`** (lines 224–347) — on activation (`Active == STAR_UNBORN`, lines 235–242),
  sets `Age = 0.0; Birthtime = All.Time;`; every subsequent call recomputes `Age = All.Time -
  Birthtime` (line 246) rather than accumulating `Age += star_timestep` incrementally. The old
  incremental scheme could drift or double-count across timestep changes or repeated `star_prep`
  calls within a step; the new one is idempotent and exact.
- **H/He feed-buffer flush** (lines 390–408, from `e906c22`, part of §1.7): at end-of-step,
  `SphP[i].GrackleSpeciesConserved(GRACKLE_SPECIES_INDEX + s) += SphP[i].StarChemFeed[s]` for
  every Grackle species, followed by `sync_primitive_from_conserved()` and zeroing the buffer.
  The metals branch immediately below was converted to the same conserved-first pattern in the
  same commit (previously primitive-first: `GasMetallicity = (...)/Mass;
  sync_conserved_from_primitive(...)`).
- `star_in()`/`star_exit()`, `star_perform_end_of_step_physics()`, `gaussian_weight(r, h)` —
  unchanged in role.

### 1.7 H and He mass injection (new since `f692626`)

A new pair of tracked species — atomic Hydrogen and Helium mass-loss (as opposed to the
aggregate `MetalsLoss`) — was threaded through the entire feedback pipeline: table storage →
interpolation → per-star feedback struct → mesh deposition → gas-cell buffer flush. Gated behind
`#if GRACKLE_CHEMISTRY >= 1`, nested inside the pre-existing `WINDS`/`SUPERNOVAE` blocks.

**Table storage (`star_tables.h`/`.c`).** New externs `HLossRate`/`HeLossRate` (3D, alongside
`MassLossRate`/`MetalsLossRate`) and `SN_HLoss`/`SN_HeLoss` (2D, alongside `SN_MassLoss`/
`SN_MetalsLoss`). `load_star_tables()` allocates, HDF5-reads (datasets `"HLossRate"`,
`"HeLossRate"`, `"SN_HLoss"`, `"SN_HeLoss"` per `(Z, M)` group), and MPI-broadcasts these exactly
like the pre-existing metals tables — **table HDF5 files must now contain these four datasets or
`load_star_tables()` fails via `terminate()`**.

**Interpolation.** `interpolate_age/mass/metallicity()` and `SN_interpolate_mass/metallicity()`
propagate `HLossRate`/`HeLossRate`/`SN_HLoss`/`SN_HeLoss` through every clamp/interpolate branch,
mirroring the existing mass/metals handling (see §1.5 for exact line numbers per function).
`star_feedback_compute()` assigns `Star.SN_HLoss`/`SN_HeLoss` (SN branch) or `Star.HLoss`/
`HeLoss = HLossRate/HeLossRate * dt` (MS branch); `units_for_feedback()` converts back to
astrophysical units (`/= All.cf_UnitMass_in_Msun`).

**Delivery (`star_feedback.c`)** — see §1.4 above for the `Feedback_Kick`/`apply_kick()`/species
index details, and the neutral-wind-vs-ionized-SN injection distinction.

**Flush (`star_update.c:390–399`)** — see §1.6 above.

**Connection to the conserved-metals refactor.** `9abee71` ("Replace primitive metallicity with
conserved metals in star formation routines") and `76f3cdb` (the analogous fix for `compute_mu`
in `src/cooling/grackle.c`) are a *separate* but related effort: swapping primitive quantities
(`GasMetallicity`, cell `Vel`) for their conserved/mass-derived counterparts (`GasMetals`,
`Momentum/Mass`) wherever a primitive field could be stale relative to the conserved state within
a hydro sub-step. `GasMetals` is `PConservedScalars[METALS_INDEX]`
(`#define GrackleSpeciesConserved(i) PConservedScalars[(i)]`, `allvars.h:1751`), with
`GRACKLE_SPECIES_INDEX = METALS_INDEX + METALS_NUMBER` (`allvars.h:361`). The H/He injection
machinery's "buffer now, merge later" pattern (`StarChemFeed[]` → flushed into
`GrackleSpeciesConserved` at end-of-step) is the same idiom `9abee71` made consistent for metals
throughout star formation — see §3 for where it also touches individual star-by-star formation.

### 1.8 SN timestepping — the "lead time" ramp (`03ab10f`)

New parameter **`All.SN_LeadTime`** (`src/main/allvars.h:1448`, registered
`src/io/parameters.c:603–604`, default `40` in `param.txt`, units Myr). In
`star_density_evaluate2()` (`star_density.c:543–563`), rather than the old scheme of blending
`SphP[i].Utherm` toward a post-SN internal energy over the fraction of a star's *entire
lifetime*, the code now computes `sn_lead_time = All.SN_LeadTime / All.UnitTime_in_Megayears`
and, only in the last `SN_LeadTime` Myr before the explosion (`time_to_sn < sn_lead_time`), ramps
a signal speed `SphP[i].Csn = max(Csn, sqrt(GAMMA*GAMMA_MINUS1*E_inject_code/Mass) * f)` where
`f = 1 - time_to_sn/sn_lead_time`. This modifies the local Courant-like signal speed used for
hydro timestep limiting, rather than mutating internal energy proxies directly (the earlier
behaviour, which was the actual "consistency" bug this commit fixes).

### 1.9 Radiative transfer (`star_radiation.h`, `star_radiation.c`, `star_radiation_tree.c`)

The wiki's description (§1.7 there) of the HEALPix ray-tracing scheme, wavebands, opacities,
H2 self-shielding, tree walk, and MPI ray exchange is still accurate and current — none of that
machinery changed in this commit range. Two things did change:

**RT timestep normalization (`a2b69a6`, `a4c7c36`).** `rt_timestep()` (`star_radiation.c:770–804`,
`#ifdef PHOTOIONIZATION`) previously looped over *all* gas particles and took the max of the
three per-species ionization rates unnormalized, with a fallback timestep of the literal `1`
(code units) when the rate was ≤ 0. It now:
- loops only over `TimeBinsHydro.ActiveParticleList` (lines 775–779) — RT timestep is computed
  only for currently hydro-active cells;
- normalizes by neutral-hydrogen fraction: `x_HI = m_HI / m_H`, where `m_H` sums the conserved
  Grackle hydrogen-bearing species (`GRACKLE_HI + GRACKLE_HII` always, plus H2/HM species under
  `GRACKLE_CHEMISTRY >= 2`, plus D species under `>= 3`), and uses `rate =
  HI_IonizationRate * x_HI` — the HeI/HeII ionization-rate terms were dropped from this
  computation entirely;
- uses a physically sensible fallback, `All.MaxSizeTimestep / All.cf_hubble_a`, instead of the
  bare `1` (line 802, the `a4c7c36` fix).

> **⚠ Known issue, not yet fixed:** the local `double rate = 0.0;` declaration that existed
> before `a2b69a6` was deleted in that commit and never reinstated. At line ~799, `rate =
> fmax(rate, SphP[i].HI_IonizationRate * x_HI);` references an undeclared identifier. This
> function will not compile under a standards-conformant C compiler when `PHOTOIONIZATION` is
> defined, until `double rate = 0.0;` is reinstated before the active-particle loop (~line 781).
> Flagging here rather than silently patching it — check with whoever owns radiation before
> fixing, in case there's a local workaround already in flight.

---

## 2. `src/blackholes/` — black hole framework

Structurally this module still closely parallels `src/stars/`: a `Bh_Particle_Data` struct
(`bh.h`) accessed via `BPP(i)`/`PPB(i)`, a `TimeBinsBh` timebin system, and the same two-phase
tree/ngb-loop pattern for density and feedback loops. The wiki's §2.1–2.7 (data model basics,
density loop, the three accretion models, swallowing, thermal/jet feedback, the imported
refinement scheme) are unchanged and still accurate. What's new since `f692626`:

### 2.1 Formation-history fields (`bh.h`, from `c0f1919`)

`Bh_Particle_Data` gained a block gated `#ifdef HALO_SEEDING` (lines 20–26):

```c
MyDouble FormationTime;         // All.Time at seeding
MyFloat  FormationMetallicity;  // donor cell's metal mass fraction at seeding; -1 if METALS is off
int      FormationChannel;      // bitmask of BH_SEED_CHANNEL_*, deliberately int (not a smaller
                                 // type) so it round-trips through the MEM_INT snapshot I/O path
MyFloat  DonorVelocity[3];      // donor cell's velocity, captured before P[ibh].Vel overwrites it
```

Written to snapshots as HDF5 fields `BFTI`/`BFMZ`/`BFCH`/`BFDV`
(`src/io/io_fields.c:841–859`, gated `#if defined(HALO_SEEDING) && defined(BLACKHOLES)`); see
[fof_halo_seeding.md](fof_halo_seeding.md) for the seeding mechanism these fields are populated
from. **Caveat**: `c0f1919`'s changes also fixed a latent bug where `BH_SEED_ON_ZERO_METALLICITY`
referenced a non-existent `SphP[].Metals` field instead of the `GasMetals` macro — see
[fof_halo_seeding.md](fof_halo_seeding.md)'s "Known limitations" section.

### 2.2 `bh_seed.c` — spawning BHs from seed events (new file, 234 lines)

Gated `#if defined(HALO_SEEDING) && defined(BLACKHOLE_SEEDING)`. Two functions:

- **`seed_black_holes_from_events(HaloSeedEvent *events, int n_events)`** (line 55) — declared in
  `fof_seeding.h`, not `bh_proto.h`, since its caller lives in `src/fof/`. Consumes a global,
  identical-on-all-tasks event array; each task filters to `events[i].DonorTask == ThisTask`,
  validates the donor gas cell hasn't gone stale since the FOF pass (index range, `Type==0`, ID
  match, not already dead), caps the seed mass (`All.BlackHoleSeedMass`) to 90% of the donor
  cell's mass with a warning if exceeded, and calls `spawn_black_hole_from_cell()`. New particle
  IDs are assigned using the **same scheme as star formation**: `calculate_maxid()` if needed,
  `MPI_Allgather` of per-task spawn counts, each task's new IDs starting at `All.MaxID + 1 +
  (sum of counts on lower-ranked tasks)`.

  Regarding `FormationChannel`: this function does **no** channel evaluation itself — it simply
  forwards `events[i].FormationChannel` (a bitmask already computed by the FOF-seeding producer
  side, per the "evaluate all channels unconditionally, don't short-circuit" logic described in
  [fof_halo_seeding.md](fof_halo_seeding.md)) straight through to the spawn function, stored
  verbatim.

- **`spawn_black_hole_from_cell()`** (line 160) — mirrors `spawn_star_from_cell()`: mass-conserving
  spawn, donor cell survives. Captures `donor_vel`/`donor_metallicity` before mutation (lines
  163–171); copies `P[igas]` into the new particle slot and overwrites `Type=5`/`Mass=seed_mass`
  (lines 173–184); scales the donor cell's `Mass`/`Energy`/`Momentum`/passive scalars by `fac =
  (Mass - seed_mass)/Mass` (lines 186–197); registers the new entry in `BhP[]`, growing it via
  `reallocate_memory_maxpartbhs()` if needed, and populates the formation-history fields (lines
  199–228, see §2.1).

### 2.3 `bh_merger.c` — BH–BH merging (new file, 259 lines, `BH_MERGER`)

**Why it exists**: nothing previously merged two black holes once their host halos combine — they
coexist indefinitely, each independently draining `SphP[].BhMassDrain` from possibly-shared gas
cells with no visibility into each other (the same root cause behind the mass-drain clamp in
§2.4). `bh_merger()` (line 58, gated `#if defined(BH_MERGER)`) closes that gap. Because the BH
count is small relative to gas particle count, it's a deliberately brute-force O(N_bh²) all-pairs
scan over a globally-gathered snapshot, executed **redundantly and identically on every MPI
task**, so no further communication is needed after the initial gather to agree on which pairs
merge.

**Algorithm:**

1. **Gather** (lines 63–108): each task builds `local_bhs[]` from its local `BhP[]` (skipping
   already-dead entries), then `MPI_Allgather`/`MPI_Allgatherv` produce a byte-identical
   `global_bhs[]` on every task. Early-out if global count < 2.
2. **Pair scoring** (lines 115–160): nested loop over all pairs, periodic-aware separation via
   `NEAREST_X/Y/Z`. Two rejection criteria:
   - **Distance**: pair rejected if `r² ≥ (All.BhMergerRadiusFactor * max(Hsml_i, Hsml_j))²`
     — the sole runtime parameter (`src/io/parameters.c:742–743`, `#ifdef BH_MERGER`).
   - **Boundedness**: pair rejected if `|Δv|² ≥ 2·G·(Mass_i+Mass_j)/r` (instantaneous
     escape-velocity test, not an orbit integration).
   Each BH tracks its own closest still-eligible partner (`best_partner[]`/`best_r2[]`).
3. **Resolve and apply** (lines 167–252): a pair merges only if it's a **mutual** nearest
   neighbour (`best_partner[j] == i` and vice versa) — since every task computes byte-identical
   arrays, every task reaches the same accept/reject decision with no further communication. The
   more massive BH survives (ties broken by lower ID); the merge is mass-summed and
   momentum-conserving (`mass_new = M_i+M_j`, `vel_new` mass-weighted mean; positions are *not*
   blended — the survivor keeps its own position). Each task applies the merge only to particles
   it owns (checked via `global_bhs[...].Task == ThisTask`), re-validating the local particle's ID
   before mutating. The loser is killed via the standard dead-particle sentinel (`Mass=0, ID=0,
   Vel={0,0,0}`) — its slot persists until the next domain decomposition; `NumMergers` (a new
   `bh.h` field, `#ifdef BH_MERGER`) on the survivor is incremented.

**LIFO bug (`fb06d35`)**: `mymalloc`/`myfree` is a strict LIFO stack allocator. An earlier version
freed `counts`/`bcounts`/`bdispls`/`local_bhs` immediately after the `Allgatherv` call, even
though `global_bhs` (allocated after all four, and needed through steps 2–3) was still on top of
the stack — this crashed with "not the last allocated block" on first call in production. Fixed
by moving those frees to the very end, in exact reverse allocation order.

**Caller contract**: collective (two collective MPI calls), every task must call it together;
intended to run once per BH-active step, *after* `bh_perform_end_of_step_physics()` has applied
that step's accretion, so mergers act on post-accretion masses. `bh_reconstruct_timebins()`
(`bh_update.c:114–118`) skips `Mass==0 && ID==0` entries so a merged-away BH drops out of
`TimeBinsBh` immediately even though its array slot persists.

### 2.4 Combined `BhMassDrain` clamp (`bh_update.c`, from `5d00766` + `e930ba3`)

`bh_perform_end_of_step_physics()` (starts line 162) gained a clamp at lines 187–207, immediately
before the pre-existing `P[i].Mass -= SphP[i].BhMassDrain;` (now line 209):

```c
if(P[i].Mass > 0)
  {
    if(SphP[i].BhMassDrain > 0.9 * P[i].Mass)
      SphP[i].BhMassDrain = 0.9 * P[i].Mass;
  }
else
  {
    if(SphP[i].BhMassDrain != 0)
      printf("BH_ACCRETION: WARNING ... corruption upstream of bh_perform_end_of_step_physics(), "
             "not from accretion draining\n", ...);
    SphP[i].BhMassDrain = 0;
  }
P[i].Mass -= SphP[i].BhMassDrain;
```

**What this fixes**: `BhMassDrain` is an accumulator that multiple BHs sharing a gas cell as an
accretion neighbour each add their own contribution to (via `bh_swallow.c`), over one step, with
no cross-visibility between BHs while accumulating — each BH previously capped only its *own*
contribution at 0.9× cell mass, so the *combined* drain from several BHs could still exceed the
cell's mass. The clamp added here is centralized and applied once per cell, after all BHs have
contributed, so the total drain is bounded regardless of how many BHs touched the cell. The
`e930ba3` follow-up fixed a sign bug in the `else` branch: naively applying the same `0.9 *
P[i].Mass` formula to an already-non-positive mass *shrinks the magnitude* of a negative value
instead of rejecting the drain — the fix instead treats `Mass ≤ 0` as evidence of upstream
corruption (not from BH accretion) and forces `BhMassDrain = 0`.

This was discovered via a real crash: a gas cell sat inside the softening length of a cluster of
5 black holes — leftover BHs from progenitor halos that had since merged — which is also the
motivating case for `bh_merger()` (§2.3); the drain clamp is the immediate fix, `BH_MERGER` is
the fix for the underlying cause.

### 2.5 Uninitialized-variable crash fixes (`bh_density.c`, `bh_swallow.c`)

`bh_density_evaluate()` (`bh_density.c:378–381`, commit `be1a178`) had `h2 = h*h;` and `hinv =
1.0/h;` commented out while the variables were still read later in the function (`if(r2 < h2)`,
`u = r*hinv`) — genuine use-of-uninitialized-value undefined behaviour, manifesting as a crash.
Fixed by restoring the two assignments.

**Not yet fully fixed**: the analogous declaration in `bh_swallow_evaluate()`
(`bh_swallow.c:222`) was uncommented in the same overall diff, but its `h2`/`hinv` *assignments*
remain commented out (lines 224–225). If `h2`/`hinv` are read later in this function the same
class of bug likely still exists here — worth a follow-up check.

### 2.6 `bh_proto.h`

New prototype `void bh_merger(void)` under `#ifdef BH_MERGER` (lines 31–34); a trailing comment
(line 55) documents that `seed_black_holes_from_events()`'s prototype deliberately lives in
`fof_seeding.h` instead. (A pre-existing, possibly-duplicate `bh_jet_density` prototype under two
different guard macros, `BH_JET_FEEDBACK` and `BLACKHOLES_FEEDBACK`, was noticed in passing —
not resolved here, flagged for whoever next touches jet feedback.)

---

## 3. Individual star-by-star formation (`src/star_formation/individual_star_formation/`)

Structure unchanged from the wiki's §3 (`individual_starbystar_formation()`, `spawn_heavy()`/
`spawn_light()`, `sf_starbystar.c`, `sfr_starbystar.c`, `sf_massdrain.c`). What changed, as part
of the same conserved-metals refactor described in §1.7 (`9abee71`):

- **`individual_star_formation.c`**: the heavy-star finalization block (lines 214–220) now
  divides `SP[i].Metallicity` by the star's final mass (`PPS(i).Mass`) once known (line 219),
  rather than that division happening inside `sf_massdrain`'s `out2particle` for only the
  imported-particle branch (removed, see below — that asymmetry meant local-mode contributions
  were never normalized there at all). `spawn_light()` (line 376) now reads
  `SP[NumStars].Metallicity = SphP[igas].GasMetals / P[igas].Mass` instead of trusting the
  (possibly stale) primitive `GasMetallicity`.
- **`sf_massdrain.c`**: `sf_massdrain_evaluate()` (lines 271–279) now accumulates velocity from
  `SphP[i].Momentum[k] / P[i].Mass` instead of the primitive `P[i].Vel[k]`, and metal mass from
  `SphP[i].GasMetals` directly instead of reconstructing `GasMetallicity × Mass` — both avoiding
  staleness relative to the conserved state within a hydro sub-step.
- **`src/star_formation/starformation.c`** — the two star-creation entry points,
  `convert_cell_into_star()` (line 355) and `spawn_star_from_cell()` (line 447), both switched
  from `SphP[...].GasMetallicity` to an on-the-fly `GasMetals / Mass` division at the exact
  instant of star formation.

The same conserved-vs-primitive idiom was applied to `compute_mu()` in `src/cooling/grackle.c`
(`76f3cdb`, outside this module but part of the same effort) since primitive values can be stale
in the hydro loop.

---

## 4. Cross-cutting infrastructure

### 4.1 BhP[]/SP[] not following P[] during FoF/SubFind output reordering

See [fof_halo_seeding.md](fof_halo_seeding.md)'s "Data-integrity fix" section for the full
writeup — this was a real correctness bug in `src/fof/fof_distribute.c` (`fof_subfind_exchange()`)
that could corrupt or crash on any snapshot write with FOF + BH/star particles active, now fixed
(`9afdd67`/`13ae729`), with a companion fix in `src/utils/allocate.c` for zeroing newly-grown
`BhP[]`/`SP[]` memory.

### 4.2 The rest, briefly (existing AREPO files extended, not new files)

- **`src/domain/*`** — MPI particle exchange and Peano-Hilbert ordering for star/BH particles;
  weights domain balancing by feedback cost.
- **`src/io/*`** — snapshot read/write fields for star/BH quantities (§2.1 above for the newest
  additions) and parameter-file options.
- **`src/mesh/voronoi/*`, `criterion_refinement.c`/`criterion_derefinement.c`** — mesh
  refinement/de-refinement triggered by star/BH proximity.
- **`src/time_integration/*`** — timestep criteria extended for star/BH activity and the
  radiative-transfer timestep (§1.9).
- **`src/cooling/*`** — the Grackle integration that H/He injection (§1.7) and SN energy
  injection ultimately feed into.

---

## 5. Known issues at HEAD (`fb06d35`)

- **`rt_timestep()` undeclared `rate`** (`star_radiation.c`, §1.9) — will not compile under
  `PHOTOIONIZATION`. Not yet fixed as of this writing.
- **`bh_swallow.c`'s `h2`/`hinv`** (§2.5) — declaration restored but assignments still commented
  out; potential latent uninitialized-read if those variables are used later in the function.
- **Duplicate `bh_jet_density` prototype** under two different guard macros in `bh_proto.h`
  (§2.6) — cosmetic, not a functional bug, but worth cleaning up.

---

## 6. Suggested reading order for a new student

1. `src/stars/star.h` and `src/blackholes/bh.h` — get the data model in your head first.
2. `src/stars/star_density.c` — see the two-phase tree-walk pattern once; every other loop in
   both modules reuses it.
3. `src/stars/star_feedback.c` (mechanical feedback) — the clearest example of the physics, with
   explicit citations (Kim & Ostriker 2015 for the cooling radius; Kimm et al. 2015 for the
   terminal momentum) to read alongside it.
4. `src/blackholes/bh_accretion.c` — pick whichever accretion model your simulations actually use
   and read only that branch.
5. `src/blackholes/bh_seed.c` + [fof_halo_seeding.md](fof_halo_seeding.md) — how BHs come into
   existence in the first place.
6. `src/blackholes/bh_merger.c` — short, self-contained, and a good example of the
   gather-everything/compute-redundantly MPI pattern used when the entity count is small.
7. `src/stars/star_radiation.c` + `star_radiation_tree.c` — leave for last; the most complex
   piece, and the one most likely to have changed again by the time you read this.

---

## Note on the wiki

The [project wiki's Technical-Reference page](https://github.com/solas-sims/arepo_solas/wiki/Technical-Reference)
predates this page and is now out of date (pinned to commit `f692626`). The plan is to replace
its content with a short pointer to this page, so there's a single source of truth that lives
alongside the code and shows up in diffs/PRs rather than drifting silently out of sync in a
separate wiki repo. That wiki edit hasn't been made yet — ask before assuming it has been.
</content>
