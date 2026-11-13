# SIDM Module for `arepo_solas`

**Status as of this writing: full v1 pipeline built, compiling, and running — density estimation validated, timestep criterion in place, scattering/kick routine implemented with cross-task partner support and a live conservation check, `DMSP[]` side-array migration complete (including restart checkpointing and snapshot I/O), and a first isolated-halo test shows the expected qualitative cusp-to-core signature. Careful quantitative validation against Vogelsberger et al. (2012) is the next milestone.**

**Update following an external review against Correa et al. (2022, TangoSIDM)**: two real gaps were found in the KDK integration hook (no position correction for the velocity change, no wake-up for a passively-kicked inactive neighbour) and fixed, gated behind `SIDM_KICK_POSITION_CORRECTION` (off by default) — see §3.10. A diagnostic mode (log detected collisions without applying them, for the isolated-halo `Gamma(r)` validation test) was added behind `SIDM_LOG_COLLISIONS` (off by default) — see §3.11. The conservation and timestep-binding diagnostics were also split/extended to directly answer two items on the TODO list below (§3.11). **All of this is unvalidated** — it compiles and links cleanly in every flag combination but has not yet been run. The kernel-overlap scattering-probability form (§3.5) was reviewed against Correa's symmetric double-kernel alternative and kept as-is; this is a closed design decision, not an open item.

This document is intended as the entry point for anyone (including future-you) picking this work back up. Read this before diving into the source. For a precise, function-by-function reference of what's implemented and where, see `IMPLEMENTATION.md` alongside this file.

---

## 1. Scope and target

- **v1 scope**: elastic-only, single-state SIDM, constant cross-section (σ/m). Validation target is the pure elastic limit of the Vogelsberger, Zavala & Loeb (2012) algorithm (as later used/extended in Vogelsberger et al. 2019, and described in detail in Valdarnini 2023, arXiv:2309.10374, which is the paper this implementation was actually built from).
- Multi-state/inelastic SIDM and a parallel FDM module (`src/fdm/`) are explicitly out of scope for v1 and deferred.
- **Codebase**: `arepo_solas`, fork of public Arepo, branch `Star_feedback_radiation`, repo at `github.com/solas-sims/arepo_solas`.
- **Working setup**: git worktree `sidm_module` (branched off `Star_feedback_radiation`) in its own directory (`git worktree add ../arepo_solas-sidm sidm_module`), so the SIDM build never shares a working tree with the baryonic-physics branch. Recommend periodic `git merge Star_feedback_radiation` into `sidm_module` to catch drift early rather than one large reconciliation later.

---

## 2. Current status — what actually works

| Piece | Status |
|---|---|
| Dedicated DM-only neighbour tree (`src/sidm/sidm_tree.c/.h`) | Built, validated (correct construction, correct multi-task behaviour, verified in both `STARS`-on and `STARS`-off builds) |
| Density/Hsml/VelDisp estimation (`sidm_density.c`) | Validated: physical density scale matches Ω_DM expectation to ~4% (cubic-spline kernel), Hsml matches independent mean-spacing prediction, converges identically on 1-task vs 4-task runs. `VelDisp` mean-subtraction bug fixed (see §3.8) |
| Timestep criterion (`timestep.c`, `get_timestep_gravity`) | Implemented, comoving-aware |
| `DMSP[]` side array | **Complete.** Struct, macros, allocation (lockstep with `P[]`), full `domain_exchange.c` mirroring of the `STARS` swap-compaction pattern, `peano.c`/`domain_rearrange.c`/`mesh/refinement.c` back-reference fixes, IC-read initial assignment, restart checkpointing. See §3.3 and §3.9 for the debugging history — several genuine gaps were found and fixed here, not just the obvious struct/macro work. |
| I/O (`io_fields.c` + `A_DMSP`) | **Complete.** Four fields (`SidmDensity`, `SidmHsml`, `SidmNumNgb`, `SidmVelDisp`) write to snapshots correctly via the new `A_DMSP` array-source dispatch (mirrors `A_S` throughout `io.c`/`read_ic.c`) |
| Scattering/kick routine (`sidm_scatter.c`) | **Local and cross-task partners both implemented.** Momentum/energy conservation verified live on real runs (both the cosmological box and the isolated halo), errors at floating-point noise level. Cross-task path compiles/links/launches cleanly but has not yet been exercised on a real multi-task run long enough to trigger actual cross-task scatters — check the `cross_task=N` count in the `SIDM_SCATTER:` diagnostic line to confirm (the line now also splits max conservation error by local vs. cross-task specifically, see §3.11). |
| KDK position-correction + neighbour wake-up (`SIDM_KICK_POSITION_CORRECTION`) | **Implemented, gated, off by default. Unvalidated.** Closes two gaps found reviewing against Correa et al. (2022): kicked positions were never corrected for the velocity discontinuity, and a passively-kicked inactive neighbour had no mechanism to be brought back onto an active timebin. See §3.10. Needs an A/B run against the flag-off behaviour before being trusted/turned on by default. |
| Collision-logging diagnostic mode (`SIDM_LOG_COLLISIONS`) | **Implemented, gated, off by default.** Logs every detected/accepted collision (both particles' positions, `r`, `p_ij`, timebin, dt) to a per-task `sidm_collisions_<task>.txt` instead of applying it, for the isolated-halo `Gamma(r)` validation test. See §3.11. The Hernquist IC generator and the `Gamma(r)` analysis script that consume this log are **not yet built** — this is the logging hook only. |
| SIDM timestep-criterion binding diagnostic (`SidmTimestepChecks`/`SidmTimestepBinding`) | **Implemented, ungated, cumulative.** Answers "does the SIDM term in `get_timestep_gravity()` ever actually bind" directly from the `SIDM_SCATTER:` line — see §3.11. Not yet read from an actual run. |
| First isolated-halo test (NFW, M=10¹² M☉/h, c=10, σ/m≈1 cm²/g) | Qualitatively healthy: monotonic core suppression over 0→120→370 Myr, outskirts (>10-20 kpc) essentially unchanged across all three snapshots (correctly localized effect), rough relaxation-timescale estimate matches the observed onset. **Not yet a quantitative match to a specific published number** — see TODO. |
| Quantitative validation against Vogelsberger's test suite | **Not started.** This is the actual "is v1 done" criterion — everything above is necessary but not sufficient. |

---

## 3. Architecture — key decisions and *why*, not just what

These are documented because each one overturned an earlier, reasonable-seeming assumption, or because a bug at that layer took real effort to track down. Skipping the "why" risks re-litigating settled questions or, worse, re-introducing a bug that was already found and fixed.

### 3.1 Dedicated tree, not gravity-tree reuse

Original plan (per Vogelsberger's own paper) was to reuse the existing gravity tree for the SIDM neighbour search, avoiding the cost of a second tree build. **This does not work in this codebase.** Under `HIERARCHICAL_GRAVITY`, `force_treebuild()` only inserts `TimeBinsGravity.NActiveParticles` particles as leaves — not the full DM population — because most timesteps only a subset of particles are active. A neighbour-density estimate needs a genuinely complete tree; the gravity tree usually isn't one. This caused non-convergent Hsml iteration in production runs before it was diagnosed and led to building `src/sidm/sidm_tree.c` as an independent, DM-only tree, structurally modelled on `src/ngbtree/ngbtree.c` (used for gas neighbour search) but:
- filtered to `Type==1` only,
- with mesh-vertex/hydro/RT-specific node fields stripped out (not applicable to static DM particles),
- rebuilt on its own cadence rather than reusing gravity's.

### 3.2 Tree rebuild cadence is tied to domain decomposition, not to a per-call basis

Second false start: initially hooked `sidm_density()` into `calculate_non_standard_physics_with_valid_gravity_tree()`, gated on `FLAG_FULL_TREE`. Two things were wrong with this, discovered in sequence:
1. `FLAG_FULL_TREE` does **not** mean "all particles active" (a doc-comment misreading, corrected after real crashes exposed it) — it fires at a coarser cadence than assumed, and irregularly relative to particle population changes.
2. More seriously: our tree-build's own domain-consistency assertion (ported from `ngb_treebuild_construct`, which checks that a particle's current position maps to a top-level domain leaf owned by the current task) is **only valid immediately after a fresh domain decomposition** — that's the only context `ngb_treebuild()` itself is ever called in. Calling `sidm_treebuild()` at any other cadence means particles can have legitimately drifted across nominal task boundaries since the last decomposition (completely normal — every other part of the codebase tolerates this), which crashed our tree build with a `STOP! ID=... should be on task=...` error.

**Resolution**: `sidm_treebuild()` and `sidm_density()` are now called explicitly, together, directly from `run.c`'s two domain-decomposition-gated blocks (same `All.HighestActiveTimeBin >= All.SmallestTimeBinWithDomainDecomposition` condition domain decomposition itself uses) — not via any gravity-tree-related hook. `sidm_treeallocate()`/`sidm_treefree()` are hooked at the same two call sites, for the same reason (`SidmTree_DomainNodeIndex` is sized off `NTopleaves`, which changes whenever decomposition changes).

**Known consequence, accepted for v1**: SIDM density only updates as often as domain decomposition runs (tunable via `All.SmallestTimeBinWithDomainDecomposition`), not every timestep. Tying rebuild-cadence to walk-cadence like this was a deliberate choice to avoid reintroducing node-level drift tracking (`ngbtree`'s `drift_node()`/`vertex_vmin`/`vertex_vmax` mechanism), which we don't have and would need if the tree were ever walked without being freshly rebuilt first.

### 3.3 `DMSP[]` side array — implemented, and harder than it looked

Originally deferred (SIDM fields lived directly on `P[]`), then implemented once "basic working version" was reached. The architectural motivation: `SP[]` (stars' own side array) is cleaner for **mixed baryon+DM production runs** since it avoids every gas/star/BH particle carrying fields it never uses — costs nothing and saves nothing in a pure-N-body test (100% DM, no non-DM population to exclude), but the real target is mixed runs.

**What this actually took, beyond the struct + macros:**
- `DM_Particle_Data` struct + `DMSP[]`/`NumDM` + `DMPS`/`PDMS` access macros (`sidm.h`), mirroring `Star_Particle_Data`/`SP[]`/`SPP`/`PPS` — but named more clearly (`PIndex` instead of reusing the confusing `PID`-as-an-array-index convention stars use).
- `allocate.c` — initial allocation + lockstep reallocation with `P[]` (no `MaxPartDM` exists, so this rides on `All.MaxPart` directly, same choice as `SidmTree_MaxPart`).
- **Full `domain_exchange.c` mirroring** of the `STARS` pattern: dedicated export buffer, per-particle export packing, local swap-removal, the actual `Sendrecv`/`Isend`/`Irecv`/`Alltoallv` calls, cross-task index reindexing. ~20 distinct touch points in one function.
- **Three additional files found only by deliberately sweeping the whole codebase** for the back-reference pattern, not files anticipated in advance: `peano.c` (Peano-key cycle-sort — runs on essentially every domain decomposition), `domain_rearrange.c` (cell elimination), `mesh/refinement.c` (gas cell refinement — potentially the *most* frequent of the three in a real AGORA run). All three had the same latent gap: shuffling the last particle into a freed slot without checking if it was DM.
- **`read_ic.c`** — the piece that caused a real bus-error crash before being found: nothing established `SIDMID`/`PIndex` for the DM population the *first* time it's read from an IC file (every case of a particle *moving* was handled; the very first assignment was not). Fixed by mirroring `STARS`'/`BLACKHOLES`' own IC-read initialization exactly (`NumDM = 0` reset, per-type counting in *both* of `read_ic.c`'s per-type reading loops, and the actual index-assignment loop).
- **`restart.c`** — checkpointing `DMSP[]`, mirroring `SP[]`'s `byten(&SP[0], NumStars*sizeof(...), modus)` pattern exactly.
- **`A_DMSP` I/O wiring** — `A_S` (the mechanism that lets `init_field()` register a field living on a side array) is not generic; it's a dedicated enum value with its own dispatch hardcoded in three separate places (`io.c`'s `init_field()` offset computation, `io.c`'s `fill_write_buffer()` write-path block, `read_ic.c`'s `empty_read_buffer()` read-path case). Needed a full `A_DMSP` mirror of all three. `get_particles_in_block()` needed *no* changes — it's genuinely generic, working off the file header's per-type counts and the `DM_ONLY` typelist bitmask rather than anything array-specific.

**A blind spot found and fixed during this work, worth knowing about generically**: my own sandbox's `Config.sh` always has `STARS` enabled, which *masked* an accidental `#ifdef SIDM` nested inside `#ifdef STARS` (in `domain_vars.c`) — the mistake only manifests in a `STARS`-off build like the isolated-halo test environment. Caught via a real build failure, then fixed by rebuilding with `STARS`/`AGORA_SF`/`USE_SFR` explicitly disabled to surface every instance of the same mistake at once (there was only the one). **Any future change touching a file that also has `#ifdef STARS` blocks should be verified in both configurations, not just the default sandbox one.**

### 3.4 Cross-section units — get this right before trusting any run

`σ/m` (cross-section per unit mass) must be supplied in **code units**, not cgs, and the conversion factor depends entirely on the active unit system:

```
[σ/m]_code = (σ/m)_cgs × UnitMass_in_g / UnitLength_in_cm²
```

**Cosmological box (Mpc/h + 10¹⁰ M☉/h convention)**: `[σ/m]_code = (σ/m)_cgs × 2.089×10⁻⁶`. An earlier test run used `SidmCrossSection = 0.1` in this convention — **≈47,900 cm²/g**, six orders of magnitude beyond any physically viable value, and the direct cause of an unexpectedly large number of timesteps (the timestep criterion was correctly reacting to the enormous effective cross-section). Corrected value: `1e-5` (≈4.8 cm²/g).

**Isolated halo (kpc + 10¹⁰ M☉ convention, no `h`)**: `[σ/m]_code = (σ/m)_cgs × 2.089`. The fiducial `1 cm²/g` benchmark (Vogelsberger et al. 2012's own primary choice, standard across the literature for isolated-halo core-formation tests) corresponds to `2.089` in this convention.

**This conversion factor is not a fixed constant of the code** — it's derived from whatever `UnitLength_in_cm`/`UnitMass_in_g` are active in `param.txt`. Recompute it explicitly any time those change.

### 3.5 Scattering algorithm and the lower-ID-initiates convention

Per-pair scattering probability (Vogelsberger et al. 2012, via Valdarnini 2023 Eq. 19):

```
P_ij = m_i · W(r_ij, h_i) · (σ/m) · v_ij · Δt_i
```

The source algorithm has each particle in a pair independently evaluate the interaction and divide by 2 (`P_i = Σ P_ij / 2`) to correct for double-counting. **We do not do this.** Instead: only the lower-ID particle in any pair evaluates and (if triggered) applies the scatter; the higher-ID particle never re-evaluates that specific pair. This structurally prevents double-evaluation in the first place, so there's nothing left for a `/2` to correct — **the `/2` is dropped, deliberately, as a derived consequence of this implementation choice**, not a verbatim transcription of the source formula. The kick itself (isotropic, CM-frame, Eq. 20) conserves momentum and kinetic energy exactly by construction — verified live on every real scatter event via a running diagnostic rather than a separate isolated test, and confirmed at floating-point noise level (~1e-15–1e-16) on both the cosmological box and the isolated halo.

**One likely typo in the secondary source, corrected**: Valdarnini (2023) states the accept condition as "P_i ≤ x" (probability less than the random draw triggers a scatter) — backwards from every standard Monte Carlo convention. Treated as a transcription error; implemented as the standard `x < P_i`.

**Reviewed against Correa et al. (2022, TangoSIDM) and kept as-is**: Correa's `g_ij` is a symmetric convolution integral of *both* particles' kernels (their Eq. 5-6, closed-form for the cubic spline in their Appendix A1.1), whereas `W(r_ij, h_i)` above uses only the initiator's own smoothing length — kernel-weighted (distance matters, not a fixed-radius uniform probability) but asymmetric, not a true overlap integral. This is a legitimate middle ground with its own citation (Valdarnini 2023 Eq. 19, the actual source used here), not an oversight. Explicitly discussed and closed as a design decision, not carried forward as an open item.

### 3.6 Cross-task partners — implemented

Originally a v1 restriction (local partners only), later implemented once the basic pipeline was validated. The core difficulty: density only ever needed to *gather* a scalar back to wherever a query originated; scattering needs to *act* on a remote particle, which the existing comm pattern was never built for.

**Implementation**: the tree walk (`sidm_scatter_walk()` in `sidm_scatter.c`) now genuinely follows pseudo-particle branches via `generic_comm_pattern`, mirroring `sidm_density.c`'s structure. A remote task, when serving a query, self-reports its own identity (`DataOut.remote_task = ThisTask`) — this can't be inferred from the comm framework's own bookkeeping. Local and remote candidates are unified into one per-particle list before accept/reject + partner selection runs. If the chosen partner is remote, both new velocities are computed on the originating task (using the remote candidate's velocity, already fetched during gathering) — the local half is applied immediately, and the remote half is delivered via a batched, one-way "set this velocity" exchange per destination task (a plain `Sendrecv` loop; not a full round-trip, since nothing needs to come back).

**Accepted approximation, discussed explicitly and agreed before implementing, not discovered after the fact**: between candidate gathering and kick delivery (both within the same `sidm_scatter()` call), a remote particle could in principle already have been kicked by a *different*, purely-local scatter on its own task. A fully rigorous fix (atomic remote claims + retry) was judged not worth the complexity for v1 — this is accepted as timestep-bounded, in the same spirit as the pre-existing (and never flagged as a bug) possibility that a purely local particle gets selected by two different lower-ID initiators in a single pass. Both are rare and bounded by the SIDM timestep criterion keeping per-step scattering probability small.

**Relies on P[] array stability within a single `sidm_scatter()` call**: a remote `P[]` index sampled during gathering must still be valid when the kick-delivery poke arrives. This holds because kicks only change velocities, never array position, and nothing reorders `P[]` mid-call — it would *not* hold if domain decomposition ran between gathering and delivery, which is why both happen back-to-back in one function.

### 3.7 Δt for the scattering probability must go through cosmological time, not a raw timebin multiply

`All.Timebase_interval`-based intervals are linear in `ln(a)` for a comoving run, **not** in physical time — `(tend-tstart)*Timebase_interval` is `Δ(ln a)`, not a duration. `sidm_dt_code_units()` in `sidm_scatter.c` goes through `get_time_difference_in_Gyr(a0, a1)` (proper cosmological-expansion-aware conversion) and then converts Gyr → code time units via `SEC_PER_YEAR`/`All.UnitTime_in_s`. Easy to get silently wrong; worth re-reading that function's comment block if anything about scattering rates looks off in a cosmological run.

### 3.8 `VelDisp` mean-subtraction fix

Originally computed as the naive `sqrt(⟨v²⟩/3)` — no mean subtraction, so any coherent/bulk flow got counted as if it were local velocity dispersion. Fixed by shipping raw first/second velocity moments (additive across local+cross-task contributions, unlike a `sqrt()`'d value) through the density pipeline and computing `var_3d = ⟨v²⟩ − ⟨v⟩²` once, after all contributions for a converged Hsml pass are in.

**A second, independent bug was found and fixed in the same pass**: the old code accumulated cross-task contributions to `VelDisp` via `P[i].SidmVelDisp += out->VelDisp` — summing already-`sqrt()`'d partial values together, which is not mathematically valid (you cannot add partial standard deviations to get a combined one). This specifically corrupted `VelDisp` for particles near a task boundary. Both bugs shared the same root cause (computing `VelDisp` prematurely, before all contributions were in) and were fixed together.

Confirmed after the fix: on a real halo-core check, `VelDisp` correctly increases with `Density` (the physically right direction) across genuinely paired samples — an earlier apparent anti-correlation turned out to be an artifact of comparing unpaired random `h5dump` draws, not a real problem.

### 3.9 Two `#ifdef` nesting bugs, both worth knowing about as a class

Two separate instances of the same mistake — inserting a new `#ifdef SIDM` block immediately after an existing `#ifdef STARS`-guarded line without checking whether that landed *inside* the `STARS` guard rather than after it:
1. `domain.h` — caught and fixed during initial review, before being handed over.
2. `domain_vars.c` — missed initially (masked by the sandbox's default `STARS`-on config, see §3.3), found via a real build failure on a `STARS`-off build, then fixed. The fix was verified by rebuilding with `STARS` off specifically to surface any *other* instances of the same mistake across every file touched that session — there was only the one remaining.

**Practical lesson, not just a fixed bug**: verify any change touching a file with adjacent `#ifdef` guards in *both* configurations you actually use, not just whichever one happens to be the default build environment.

### 3.10 KDK hook gaps closed: post-kick position correction + neighbour wake-up, gated behind `SIDM_KICK_POSITION_CORRECTION`

Found reviewing the KDK integration hook against Correa et al. (2022, Sec 2.3)'s SWIFT implementation. Their scheme interleaves the SIDM kick *between* the drift and the final half-kick of a KDK step, and corrects positions with an explicit backward/forward `D(dt/2)` pair straddling that point; they also explicitly warn that an active particle kicking an otherwise-inactive neighbour must wake that neighbour for the following timestep.

**Our hook is structurally different**: `sidm_scatter()` is called from `do_gravity_hydro.c` *after* this timebin's gravity kick has already fully completed (not between a split half-kick), so there's no clean midpoint to straddle the way Correa's scheme does. Before this fix, the kick functions (`sidm_apply_kick_local`/`sidm_apply_kick_cross_task`, and the cross-task kick-delivery receive loop) only ever touched `P[].Vel` — positions were never corrected for the velocity discontinuity, and a passively-kicked particle's `TimeBinGrav` was never touched, so it could carry a stale, too-long timestep for many steps after being kicked.

**Fix, both gated behind `SIDM_KICK_POSITION_CORRECTION` (off by default)**:
- `sidm_apply_position_correction()` retroactively corrects `P[].Pos` by `(v_new - v_old) * dt_drift_half`, using `get_drift_factor()` (the same comoving-aware mechanism `drift_particle()` itself uses, not a raw `v*dt/2`) over half of the *initiator's own* timebin step. This is an **adaptation** of Correa's scheme to this codebase's different hook point, not a transcription — the choice of "half the initiator's step" as the correction window is a judgment call, not a derived quantity.
- `sidm_wake_particle()` forces a kicked-but-otherwise-inactive particle's `TimeBinGrav` down to `All.LowestActiveTimeBin` via `timebin_move_particle()`, so it's picked up by the ordinary gravity/timestep cascade at the next global step. Deliberately only touches the particle's persistent bin membership, not the in-progress `TimeBinsGravity` active-list snapshot (unsafe to mutate mid-iteration).
- For the cross-task path, a `dv[3]` field was added to `sidm_kick_poke` so the *receiving* task (which owns the kicked particle) can apply both fixes itself — the originating task can't do either for a particle it doesn't own.

**Off by default, deliberately**: this changes physics (positions, and which timebin a particle sits on), so the exact pre-fix behaviour needs to stay available for an A/B comparison against this fix on the isolated-halo test before it's trusted. See TODO.

### 3.11 Validation diagnostics: conservation split by local/cross-task, timestep-binding counter, and collision-only logging mode

Three additions, all aimed at making the existing TODO items ("confirm cross-task scattering is exercised and conserves cleanly", "verify the timestep criterion actually binds") answerable from a real run's output rather than inferred indirectly.

**Conservation error, split**: `sidm_check_conservation()` now takes an `is_cross_task` flag and tracks momentum/energy error maxima *separately* for local vs. cross-task kicks (`sidm_max_momentum_error_local`/`_cross_task` etc.), rather than one pooled max. A pooled max can't tell you whether cross-task kicks conserve cleanly if they're a small minority of a run's total scatters — the previous single-max version was silently unable to answer that TODO item even with real cross-task data in hand. Printed as separate columns in the `SIDM_SCATTER:` line, with a note appended when `cross_task=0` so the columns aren't misread as "cross-task conserves perfectly" when really no cross-task event has happened yet.

**Timestep-binding counter**: new cumulative (not per-step) globals `SidmTimestepChecks`/`SidmTimestepBinding` (`allvars.h`/`.c`), incremented in `get_timestep_gravity()`'s SIDM block (`timestep.c`) every time the term is evaluated / actually turns out to be the smallest. Reduced and printed in the `SIDM_SCATTER:` line as `timestep_binding=N/M (X%)` — printed from `sidm_scatter.c`, not `timestep.c` itself, since the latter runs once per particle per timestep reassignment, far too often for its own diagnostic print. Ungated (a few integer increments, no I/O) since it answers an existing validation question rather than changing physics.

**Collision-only logging mode, `SIDM_LOG_COLLISIONS` (off by default)**: for the isolated-halo `Gamma(r)` validation test (§6 TODO, and the eventual analysis script this feeds). When set, `sidm_scatter()`'s dispatch calls `sidm_log_collision()` **instead of** the normal kick path — a genuinely separate branch, not a "kick then undo" — so `P[]` is provably untouched and the halo's density profile stays fixed and known for the whole run. Logs both particles' actual positions (not just the initiator's): for a local candidate, the partner's position is read live from `P[cand->index].Pos`; for a remote candidate, it's shipped through the existing candidate-gathering `data_out` response (a new `pos[MAX_REMOTE_RESPONSE][3]` field, mirroring how `vel` was already shipped), gated the same way so the normal physics path pays nothing extra for it. Output goes to a per-task `sidm_collisions_<task>.txt` (mirrors the pre-existing `FdDetailed`/`DETAILEDTIMINGS` per-task-file pattern in `logs.c`, since collisions are detected wherever the initiating particle lives, not just on the root task).

**Still needed for the actual validation test** (not built by this diagnostic hook alone): the Hernquist-only isolated-halo IC (reuse `examples/isolated_galaxy_collisionless_3d`'s GalIC setup with `MD=MB=N_DISK=N_BULGE=0` — not yet confirmed this degenerates cleanly to halo-only), and the Python script that bins the logged collisions radially and compares against the analytic `Gamma(r)` prediction.

---

## 4. File manifest

See `IMPLEMENTATION.md` for a complete, function-by-function reference. Summary:

```
src/sidm/
  sidm.h            — module header; DM_Particle_Data struct, DMSP[]/NumDM, DMPS/PDMS macros,
                      SIDM_TIMESTEP_SAFETY_FACTOR
  sidm_tree.h/.c    — dedicated DM-only tree: allocation, construction, domain-topleaf integration
  sidm_density.c    — Hsml-iteration density/velocity-dispersion estimator, cubic-spline kernel,
                      mean-subtracted VelDisp
  sidm_scatter.c    — Monte Carlo pairwise scattering/kick routine, local + cross-task partners

Touched outside src/sidm/:
  src/main/allvars.h         — All. fields, DM_ONLY, A_DMSP, IO_SIDM_* tags, particle_data.SIDMID
  src/main/run.c             — sidm_treeallocate/free + sidm_treebuild + sidm_density call sites
  src/io/io.c                — A_DMSP dispatch (init_field offset, fill_write_buffer write path)
  src/io/io_fields.c         — 4 field registrations (SidmDensity/Hsml/NumNgb/VelDisp), now A_DMSP-based
  src/io/read_ic.c           — A_DMSP read-path dispatch; DM SIDMID/PIndex initial assignment from IC
  src/io/restart.c           — DMSP[] checkpointing
  src/io/parameters.c        — SidmDesNumNgb/SidmDesNumNgbDev/SidmCrossSection param.txt registration
  src/utils/allocate.c       — DMSP[] initial allocation + lockstep reallocation with P[]
  src/utils/tags.h           — TAG_DMDATA
  src/domain/domain.h        — toGoDM/toGetDM declarations
  src/domain/domain_vars.c   — toGoDM/toGetDM allocation/free
  src/domain/domain_counttogo.c — DM per-particle counting pass
  src/domain/domain_exchange.c  — full DMSP[] cross-task migration (buffer, packing, swap-removal,
                                  Sendrecv/Isend/Irecv/Alltoallv, reindexing)
  src/domain/domain_rearrange.c — DM back-reference fixup (cell elimination)
  src/domain/peano.c            — DM back-reference fixup (Peano-key cycle-sort)
  src/mesh/refinement.c         — DM back-reference fixup (gas cell refinement)
  src/time_integration/timestep.c        — SIDM timestep criterion in get_timestep_gravity;
                                            SidmTimestepChecks/SidmTimestepBinding increments (§3.11)
  src/time_integration/do_gravity_hydro.c — sidm_scatter() call site, after gravity's per-timebin kick
  src/io/logs.c              — FdSidmCollisions open/close (SIDM_LOG_COLLISIONS, §3.11)
  Config.sh                  — SIDM flag; SIDM_KICK_POSITION_CORRECTION (§3.10) and
                                SIDM_LOG_COLLISIONS (§3.11), both off by default
  Makefile                   — sidm/*.o build rules
```

---

## 5. Required reading

Before making further changes, read (in this order):

1. **This file**, then `IMPLEMENTATION.md` for the precise code reference.
2. `src/sidm/sidm_tree.h` header comment — why this tree exists separately from both the gravity tree and `ngbtree`.
3. `src/sidm/sidm_density.c` header comment + driver comment — domain-decomposition cadence tie-in, mean-subtraction fix.
4. `src/sidm/sidm_scatter.c` header comment — lower-ID-initiates derivation, cross-task partner mechanism, the accepted staleness approximation.
5. `src/sidm/sidm.h` — struct/macro naming, and why `PIndex` was chosen over reusing `SP[]`'s confusing `PID`-as-array-index convention.
6. Valdarnini (2023), arXiv:2309.10374, Section 2.3 — the actual source of the scattering algorithm implemented here. Not the primary Vogelsberger et al. (2012) paper itself — worth going to the primary source if exact fidelity ever matters.
7. Vogelsberger, Zavala & Loeb (2012), MNRAS 423, 3740 — original elastic algorithm (not yet independently verified against; see above).

---

## 6. TODO / known limitations, roughly in priority order

### Blocking "v1 is actually done"
- [ ] **Careful quantitative validation against Vogelsberger et al. (2012)'s own test suite.** The isolated-halo run so far is a genuine, encouraging qualitative match (monotonic core suppression, correctly localized to small radii, rough relaxation-timescale agreement) — but "looks right" and "quantitatively matches a specific published number at a specific time" are different bars. This is the actual "done" criterion. The `SIDM_LOG_COLLISIONS` diagnostic mode (§3.11) is a step toward this but is not itself the validation — still need the Hernquist-only IC and the `Gamma(r)` analysis script.
- [ ] Confirm cross-task scattering is genuinely exercised on a real run (check `cross_task=N` in the diagnostic) and that conservation still holds cleanly for cross-task kicks specifically, not just local ones. **The diagnostic can now answer this directly** (`SIDM_SCATTER:` line splits `max_p/ke_rel_err` by local/cross-task, §3.11) — just needs an actual ≥2-task run to read the output from.
- [ ] Verify the timestep criterion actually *binds* somewhere in a real run. **The diagnostic can now answer this directly** (`timestep_binding=N/M` in the `SIDM_SCATTER:` line, §3.11) — same as above, needs a real run.
- [ ] **New**: A/B the `SIDM_KICK_POSITION_CORRECTION` fix (§3.10) against the pre-fix behaviour on the isolated-halo test. The fix compiles and links cleanly in every flag combination but has not been run — don't assume it's a net improvement (or even harmless) without actually comparing density profiles and the conservation diagnostic flag-on vs. flag-off.
- [x] ~~Kernel-overlap scattering-probability form (external review item)~~ — reviewed against Correa et al. (2022)'s symmetric double-kernel alternative and kept as-is; closed as a design decision, see §3.5.

### Architecture, deferred by deliberate choice (not urgent, but real)
- [ ] Non-`HIERARCHICAL_GRAVITY` build path — `sidm_scatter()` is not wired into `do_gravity_hydro.c`'s non-hierarchical branch (structurally different active-list handling; not exercised by any current build).
- [ ] A fully rigorous cross-task staleness fix (atomic remote claims + retry) — explicitly not done for v1, see §3.6.

### Smaller, worth tracking
- [ ] `SIDM_TIMESTEP_SAFETY_FACTOR` is hardcoded (`0.1`) rather than a `param.txt` entry — **and, separately, its actual value is an open question**: Vogelsberger/Correa's own `kappa` is `1e-2`, a full order of magnitude tighter than the `0.1` used here. Not yet decided whether to tighten it or justify `0.1` explicitly; needs testing before either call.
- [ ] `MAX_CAND_TOTAL` (128) / `MAX_REMOTE_RESPONSE` (16) in `sidm_scatter.c` are generous fixed caps, not dynamically sized — silently drop excess candidates rather than erroring if ever exceeded (worth an occasional check that this isn't happening in practice).
- [ ] The pre-existing `TAG_BHDATA`/`TAG_STARDATA` mismatch found in the *original* `STARS` code (`domain_exchange.c`'s `Sendrecv` for star exchange uses different send/recv tags — a real latent bug, not ours, not fixed by us) — flagged to Nick (owns that code) separately.
- [ ] The `Config_AGORA.sh` Grackle discrepancy noted very early on (the version pasted at project start had `USE_GRACKLE`/`METALS` disabled; the version on GitHub's `Star_feedback_radiation` branch has them active) was flagged once and never resolved.
- [x] ~~`domain_exchange.c`'s local swap-removal logic has a generic `else` fallback branch (for whichever particle ends up shuffled into a freed slot) that never got an SIDM/`DMPS` cross-check, unlike the GAS/STARS/BLACKHOLES branches which each explicitly fix up their own type's back-reference there.~~ Fixed: added the missing `#ifdef SIDM if(P[NumPart-1].Type == 1) DMPS(NumPart-1).PIndex = n; #endif` block, matching the STARS/BLACKHOLES pattern immediately above it in the same branch. Compiles/links clean; not yet exercised on a real multi-task run (this branch only fires for a particle type with no dedicated case of its own, e.g. gas after the Type==0 branch above already handled it — worth confirming it's actually reachable for a Type==1 particle in practice, not just structurally correct).

---

## 7. Practical notes for whoever runs this next

- Two validation environments now exist: the original 5 Mpc/h, 64³-particle DM-only cosmological box (cheap, useful for infrastructure bugs, says nothing about physics correctness) and a proper isolated NFW halo (M=10¹² M☉/h, c=10, 10⁵ particles within R_vir, non-cosmological) — the latter is the one that actually matters for physics validation.
- Every multi-task bug found so far (the `out2particle` cross-task corruption bug, the domain-mismatch tree crash, the `#ifdef` nesting bugs) was invisible on single-task or default-config runs. **Always sanity-check anything new on ≥2 tasks, and in both `STARS`-on/off configurations, before trusting it.**
- If density/Hsml/NumNgb ever look wrong again, the diagnostic pattern that worked well throughout this project: (1) check aggregate MPI-reduced min/max/mean first, (2) if something looks off, track a **fixed set of particle IDs across consecutive iterations** rather than fresh samples each time (a fresh sample can accidentally show only "easy" cases and hide the real problem), (3) once a misbehaving particle is identified, trace local-only vs. cross-task contributions separately.
- When comparing values across snapshots by hand (e.g. via `h5dump`), **make sure you're comparing the same particles**, not independent random draws — an apparent `Density`/`VelDisp` anti-correlation early on turned out to be exactly this mistake, not a real physics problem.
