# SIDM Module for `arepo_solas`

**Status as of this writing: basic end-to-end pipeline built and compiling — density estimation validated, timestep criterion in place, scattering/kick routine implemented with a live conservation check, quantitative validation against the literature not yet done.**

This document is intended as the entry point for anyone (including future-you) picking this work back up. Read this before diving into the source.

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
| Dedicated DM-only neighbour tree (`src/sidm/sidm_tree.c/.h`) | Built, compiling, validated (correct construction, correct multi-task behaviour) |
| Density/Hsml/VelDisp estimation (`sidm_density.c`) | Validated: physical density scale matches Ω_DM expectation to ~4% (cubic-spline kernel), Hsml matches independent mean-spacing prediction, converges identically on 1-task vs 4-task runs |
| Timestep criterion (`timestep.c`, `get_timestep_gravity`) | Implemented, comoving-aware, not yet empirically checked for whether it actually binds in a real run |
| I/O (`io_fields.c`) | Four fields (`SidmDensity`, `SidmHsml`, `SidmNumNgb`, `SidmVelDisp`) write to snapshots correctly |
| Scattering/kick routine (`sidm_scatter.c`) | Implemented, compiling, running; momentum/energy conservation diagnostic added but **results not yet reviewed** |
| Quantitative validation against Vogelsberger's test suite | **Not started.** This is the actual "is v1 done" criterion — everything above is necessary but not sufficient. |

---

## 3. Architecture — key decisions and *why*, not just what

These are documented because each one overturned an earlier, reasonable-seeming assumption. Skipping the "why" risks re-litigating settled questions or, worse, re-introducing a bug that was already found and fixed.

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

### 3.3 `P[]` fields now, `DMSP[]` side array later (deliberately deferred)

SIDM per-particle fields (`SidmDensity`, `SidmHsml`, `SidmVelDisp`, `SidmNumNgb`, `SidmLastScatterTime`, `SidmScatterFlag`) currently live directly on the shared `particle_data` struct (`P[]`), the same array every gas cell, star, and BH also uses.

**This was a deliberate, discussed trade-off**, not an oversight. The alternative — a dedicated `DMSP[]` side array, mirroring how `SP[]` works for stars (`#define PPS(i) P[SP[i].PID]`) — is architecturally cleaner for **mixed baryon+DM production runs** (avoids wasting memory on gas/star/BH particles that never need these fields), but costs nothing and saves nothing in a **pure N-body test** (100% DM particles — no non-DM population to exclude, so the side array's back-reference field is pure overhead in that specific case). Given the plan is "basic working version first, polish later," `P[]`-based fields were kept for now. **`DMSP[]` is the single largest piece of deferred architecture work** — see TODO list.

**Direct consequence**: `domain_exchange.c` doesn't need a `Type==1`-specific compaction branch yet (SIDM fields migrate for free along with the rest of `P[]`), but **will** the moment `DMSP[]` exists — that side array needs its own swap-compaction logic during MPI particle migration, mirroring the `STARS`/`BLACKHOLES` branches already in `domain_exchange.c`.

### 3.4 Cross-section units — get this right before trusting any run

`σ/m` (cross-section per unit mass) must be supplied in **code units**, not cgs. For this codebase's Mpc/h + 10¹⁰ M☉/h convention:

```
[σ/m]_code = (σ/m)_cgs × 2.089×10⁻⁶
(σ/m)_cgs  = [σ/m]_code × 4.787×10⁵
```

An earlier test run used `SidmCrossSection = 0.1`, which is **≈47,900 cm²/g** — six orders of magnitude beyond any physically viable value, and the direct cause of an unexpectedly large number of timesteps (the timestep criterion was correctly reacting to the enormous effective cross-section). Corrected value used in subsequent testing: `1e-5` (≈4.8 cm²/g, comfortably inside the literature's fiducial range of ~0.1–10 cm²/g).

**If you ever change `UnitLength_in_cm`/`UnitMass_in_g` in `param.txt`, this conversion factor changes too** — it is not a fixed constant of the code, it's derived from whatever units convention is active.

### 3.5 Scattering algorithm and the lower-ID-initiates convention

Per-pair scattering probability (Vogelsberger et al. 2012, via Valdarnini 2023 Eq. 19):

```
P_ij = m_i · W(r_ij, h_i) · (σ/m) · v_ij · Δt_i
```

The source algorithm has each particle in a pair independently evaluate the interaction and divide by 2 (`P_i = Σ P_ij / 2`) to correct for double-counting. **We do not do this.** Instead: only the lower-ID particle in any pair evaluates and (if triggered) applies the scatter; the higher-ID particle never re-evaluates that specific pair. This structurally prevents double-evaluation in the first place, so there's nothing left for a `/2` to correct — **the `/2` is dropped, deliberately, as a derived consequence of this implementation choice**, not a verbatim transcription of the source formula. This matters if you're ever cross-checking implementation against the paper line-by-line.

### 3.6 Local-partners-only (v1 restriction)

A scattering partner must be a genuinely local `P[]` particle — candidates reached only via the tree's pseudo-particle/cross-task branch are skipped entirely. Applying a symmetric kick to a cross-task partner would need new "send an action to be applied on another task" MPI infrastructure that doesn't exist yet (our existing comm pattern only ever gathers a scalar result back, it doesn't modify remote particle state). This is a real, bounded approximation — particles near a task boundary have a slightly undercounted candidate pool — deferred to the same "polish later" bucket as `DMSP[]`.

### 3.7 Δt for the scattering probability must go through cosmological time, not a raw timebin multiply

`All.Timebase_interval`-based intervals are linear in `ln(a)` for a comoving run, **not** in physical time — `(tend-tstart)*Timebase_interval` is `Δ(ln a)`, not a duration. `sidm_dt_code_units()` in `sidm_scatter.c` goes through `get_time_difference_in_Gyr(a0, a1)` (proper cosmological-expansion-aware conversion) and then converts Gyr → code time units via `SEC_PER_YEAR`/`All.UnitTime_in_s`. Easy to get silently wrong; worth re-reading that function's comment block if anything about scattering rates looks off in a cosmological run.

---

## 4. File manifest

```
src/sidm/
  sidm.h            — module header (mirrors src/stars/star.h's role); SIDM_TIMESTEP_SAFETY_FACTOR constant
  sidm_tree.h/.c    — dedicated DM-only tree: allocation, construction, domain-topleaf integration
  sidm_density.c    — Hsml-iteration density/velocity-dispersion estimator, cubic-spline kernel
  sidm_scatter.c    — Monte Carlo pairwise scattering/kick routine

Touched outside src/sidm/:
  src/main/allvars.h        — particle_data fields, All. fields (SidmDesNumNgb/Dev, SidmCrossSection),
                              DM_ONLY type constant, IO_SIDM_* tags
  src/main/run.c            — sidm_treeallocate/free + sidm_treebuild + sidm_density call sites
  src/main/proto.h          — NOT touched for SIDM (module prototypes deliberately live in sidm.h,
                              matching the codebase's own convention of not duplicating them here)
  src/io/io_fields.c        — 4 field registrations (SidmDensity/Hsml/NumNgb/VelDisp)
  src/io/parameters.c       — SidmDesNumNgb/SidmDesNumNgbDev/SidmCrossSection param.txt registration
  src/time_integration/timestep.c        — SIDM timestep criterion in get_timestep_gravity
  src/time_integration/do_gravity_hydro.c — sidm_scatter() call site, after gravity's per-timebin kick
  Config.sh                 — SIDM flag
  Makefile                  — sidm/*.o build rules
```

---

## 5. Required reading

Before making further changes, read (in this order):

1. **This file.**
2. `src/sidm/sidm_tree.h` header comment — explains why this tree exists separately from both the gravity tree and `ngbtree`, and what was deliberately left out of the node struct.
3. `src/sidm/sidm_density.c` header comment + `sidm_density()` driver comment — explains the domain-decomposition cadence tie-in.
4. `src/sidm/sidm_scatter.c` header comment — explains the lower-ID-initiates derivation and the local-partners-only restriction.
5. Valdarnini (2023), arXiv:2309.10374, Section 2.3 — the actual source of the scattering algorithm implemented here. Not the primary Vogelsberger et al. (2012) paper itself — worth going to the primary source if exact fidelity ever matters (see TODO).
6. Vogelsberger, Zavala & Loeb (2012), MNRAS 423, 3740 — original elastic algorithm (not yet independently verified against; see above).

---

## 6. TODO / known limitations, roughly in priority order

### Blocking "v1 is actually done"
- [ ] **Review the momentum/energy conservation diagnostic output** from the current running test — added but results not yet confirmed. Expect errors at floating-point noise level (~1e-14); anything above ~1e-10 is a real bug.
- [ ] **Quantitative validation against Vogelsberger's elastic-limit test suite** — an isolated halo core-formation test, not the cosmological box used for development so far. This is the actual "done" criterion for v1; nothing above substitutes for it.
- [ ] Verify the timestep criterion actually *binds* somewhere in a real run (i.e. `dt_sidm` is sometimes the limiting factor) — if it never does, understand why before trusting it's doing anything.

### Architecture, deferred by deliberate choice (not urgent, but real)
- [ ] `DMSP[]` side array (Section 3.3) — the single biggest deferred piece. Needs: struct + `PID`/`SID`-style bidirectional cross-reference, `domain_exchange.c` `Type==1` swap-compaction branch (mirroring `STARS`/`BLACKHOLES`), `allocate.c` wiring, restart I/O.
- [ ] Cross-task scattering (Section 3.6) — needs new MPI infrastructure to apply a kick to a remote particle, not just gather a scalar.
- [ ] `VelDisp` mean-subtraction — still the naive `sqrt(⟨v²⟩)` formula, not `sqrt(⟨(v−v̄)²⟩)`. Confirmed not to break the qualitative density-dispersion trend (checked against real, correctly-paired data), but will bias absolute scattering rates. Fix once v1 is validated, before trusting quantitative results.
- [ ] Non-`HIERARCHICAL_GRAVITY` build path — `sidm_scatter()` is not wired into `do_gravity_hydro.c`'s non-hierarchical branch (that branch's active-list structure differs enough that writing an untested call there wasn't worth the risk; also not exercised by any current build).

### Smaller, worth tracking
- [ ] `SIDM_TIMESTEP_SAFETY_FACTOR` is hardcoded (`0.1`) rather than a `param.txt` entry — promote if tuning turns out to matter.
- [ ] `MAX_SCATTER_CANDIDATES` fixed at 512 in `sidm_scatter.c` — comfortably above the ~32-neighbour target, but silently drops excess candidates rather than erroring if ever exceeded.
- [ ] Restart-safety of new `particle_data` fields assumed automatic (generic struct serialization) but never explicitly confirmed by inspecting a restart file.
- [ ] The `Config_AGORA.sh` Grackle discrepancy noted very early on (the version pasted at project start had `USE_GRACKLE`/`METALS` disabled; the version on GitHub's `Star_feedback_radiation` branch has them active) was flagged once and never resolved — worth reconciling before any run that needs to match the "real" branch state.

---

## 7. Practical notes for whoever runs this next

- Development/testing has used a 5 Mpc/h, 64³-particle DM-only cosmological box. This is **not** the validation target — it's been useful for shaking out infrastructure bugs (tree correctness, MPI determinism, cadence issues) precisely because it's cheap to run and easy to reason about, but it says nothing about whether the *physics* is quantitatively right.
- Every multi-task bug found so far (the `out2particle` cross-task corruption bug, the domain-mismatch tree crash) was invisible on single-task runs. **Always sanity-check anything new on ≥2 tasks before trusting it.**
- If density/Hsml/NumNgb ever look wrong again, the diagnostic pattern that worked well throughout this project: (1) check aggregate MPI-reduced min/max/mean first, (2) if something looks off, track a **fixed set of particle IDs across consecutive iterations** rather than fresh samples each time (a fresh sample can accidentally show only "easy" cases and hide the real problem), (3) once a misbehaving particle is identified, trace local-only vs. cross-task contributions separately.
