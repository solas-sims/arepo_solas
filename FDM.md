# FDM (Fuzzy Dark Matter) Module for `arepo_solas`

**Status as of this writing: Phase 1 complete AND wired into Arepo's actual run loop, validated via a real soliton stability test (multi-resolution: N=64/128/256, analysis tooling built and debugged against real output). Phase 2a complete: bidirectional star↔wavefunction coupling (mesh→star potential/force interpolation, star→mesh mass deposition, both hooked into the real simulation loop and Arepo's own gravity accumulation), validated in isolation per-function AND via a real end-to-end run. No gas coupling (2b) or BH coupling (2c) yet.**

This document is the entry point for anyone (including future-you) picking this work back up. Read this before diving into the source. For a precise, function-by-function reference, see `IMPLEMENTATION.md` alongside this file.

---

## 1. Why this module exists, and what it's actually for

The broader goal is a **pluggable dark matter backend** within `arepo_solas`: CDM (ordinary N-body), SIDM (already built, see the separate SIDM module), and FDM (this module) as interchangeable DM physics under the *same* baryonic pipeline, IC conventions, and analysis tooling — enabling direct, apples-to-apples comparison between DM models on identical galaxy-formation setups, something no single external codebase currently offers in one place.

**Specific science target** (not a generic "implement FDM" exercise): how a solitonic core and granule-driven density fluctuations affect the assembly and long-term dynamical stability of a central galaxy's stellar population, plus stream/shell morphology from tidal disruption. This requires genuine wave physics (solitonic cores, interference granules) — a cheaper "quantum pressure as an effective N-body force" approximation was considered and explicitly rejected for this specific science target, since it cannot produce real interference structure or granule-driven heating (see §3.1).

**Deliberately checked for novelty before committing significant time**: searching turned up a closely related, already-published line of work (the "Galaxy Formation with wave/fuzzy dark matter" / BECDM series, using Arepo with a similar baryonic physics lineage) that appears to have already built bidirectional FDM+full-baryonic-physics coupling in Arepo. This module is being built anyway, because (a) we don't have access to that code, (b) the specific science question here (granule-driven heating/stability, streams/shells) doesn't appear to be what those papers targeted, and (c) having FDM alongside our own validated SIDM module in one codebase, under one pipeline, has standalone value regardless of whether the base coupling mechanism is novel. Worth knowing this context before assuming the whole effort is unprecedented.

---

## 2. Scope: what Phase 1 is, and the staged plan beyond it

**Phase 1 (this document's subject)**: a standalone, non-cosmological, DM-only pseudo-spectral wavefunction solver, following May & Springel (2021)'s numerical method (split-step kick-drift-kick, Eqs. 17-20 of that paper) but reusing Arepo's own existing FFT infrastructure rather than the from-scratch implementation that paper describes.

**Central-galaxy focus, not cosmological boxes**: May & Springel's own cosmological-box costs are firmly impractical (their largest run: >7×10⁶ CPU-hours). Restricting to a single central galaxy — not a full halo out to the virial radius, and definitely not a cosmological volume — changes the cost by many orders of magnitude for two compounding reasons: cost scales as `N⁵` at fixed resolution, and a central-galaxy box is far smaller than a cosmological one, so `N` itself shrinks proportionally. See §3.2 for the actual sizing exercise and a genuinely counterintuitive finding from it.

**Staged plan beyond Phase 1** (each phase a deliberate derisking step, not just a checklist):
- **Phase 2a — COMPLETE**: couple the wavefunction to *stars only* (no gas, no star formation) — stars are already collisionless particles in Arepo, structurally the easiest first coupling target, and this establishes the "particle feels a mesh-sourced potential, mesh feels particle mass back" pattern that any future coupling (including a hybrid zoom's N-body region, see below) would need generally. See §3.9-3.11 for the architecture and hard-won lessons; §6 for what's still a known simplification within this phase (sync-point-only coupling, no gas/BH).
- **Phase 2b**: add gas — the genuinely hard part, since Arepo's moving Voronoi mesh needs to both source and feel the wavefunction's potential.
- **Phase 2c**: add the central black hole — incremental on top of 2a/2b, but BH dynamics inside a solitonic core is itself a live, actively-studied question, not just a completeness item.
- **Phase 3**: actual science — stability/assembly diagnostics, stream/shell analysis.
- **Phase 5 (deferred, scoped but not started)**: a hybrid zoom capability — N-body at large scales, full wave dynamics in a nested high-resolution patch, mirroring the literature's own AMR-based approach (Schwabe/Veltmaat/Niemeyer, GAMER-2) but *without* building genuine AMR (a much larger undertaking than this module; Arepo's moving-mesh architecture has no native AMR to build on, which is likely *why* May & Springel built the expensive uniform-mesh version in the first place rather than the cheaper hybrid approach other groups use). A fixed or periodically-recentred nested patch — not dynamically self-refining — was identified as a tractable middle ground, deferred until Phase 1-3 validate the core machinery. Phase 2a's transpose-based, task-ownership-driven communication pattern (§3.10) is directly relevant groundwork here, not just Phase 2a-specific plumbing.

**Deliberately not pursued**: full dynamic AMR (too large an undertaking relative to the payoff at this stage); the "quantum pressure as effective force" approximation (rejected for the reason in §1).

---

## 3. Architecture — key decisions and *why*

### 3.1 Why full wave dynamics, not an N-body approximation

Considered explicitly and rejected: adding an extra force term to ordinary DM particles approximating the leading-order quantum-pressure effect (validated in the literature as a comparison method, e.g. Veltmaat & Niemeyer 2016). Dramatically cheaper, and would have fit more naturally alongside the existing SIDM particle infrastructure. Rejected because the specific science target (granule-driven heating, genuine interference structure) is fundamentally about real, stochastic wave interference — no force-based approximation produces that, by construction.

### 3.2 Cost sizing — a genuinely counterintuitive finding, worth remembering

Central-galaxy scoping is not a minor convenience; it changes feasibility by orders of magnitude. Worked example (`mc²=10⁻²² eV`, kpc+10¹⁰M☉+km/s units, calibrated against May & Springel's own quoted cost at their published resolution): a `L=50 kpc`, `N=256` box for a Milky-Way-mass halo (`M=10¹² M☉/h`) comes to roughly **5,600 CPU-hours** — a ~1,250× reduction from their reference run, from a combination of smaller box (`L²` term) and smaller mesh (`N⁵` term).

**The counterintuitive part**: naively, a *lower*-mass halo should be cheaper still. It isn't, necessarily. For a dwarf-scale halo (`M=10⁸ M☉/h, c=30`), the relevant resolution constraint flips: the pure velocity/phase-aliasing criterion (Eq. 22/23 of the paper) becomes far *less* stringent (lower virial velocities → longer de Broglie wavelength → coarser allowed `Δx`), but the halo's own structural scale (NFW scale radius, soliton core) is far *smaller* in absolute terms — and resolving that structural scale, not the velocity criterion, becomes the binding constraint. A first pass using the velocity criterion alone for this halo gave a nonsensical `N=4` result; correctly using `min(velocity-criterion Δx, structural Δx)` gave a *much* higher cost than the massive-halo case, until the box size was also scoped consistently (just the central region, not the whole halo) — at which point cost came back down to a comparable order of magnitude. **Lesson for any future halo choice**: check both criteria explicitly, and scope the box to the region actually being studied, not the whole halo, for every mass regime — don't assume the massive-halo intuition transfers to lower masses.

### 3.3 Reusing Arepo's existing FFT infrastructure — what's genuinely reusable, what isn't

Arepo's own TreePM long-range gravity solver (`src/gravity/pm/`) already has a distributed, MPI-parallel FFT engine, including a genuine complex-to-complex path (`my_slab_based_fft_c2c`, `my_slab_transpose`) that existed in the codebase but was **completely unused elsewhere** before this module. This is a substantial, already-correct piece of infrastructure (the distributed data layout and MPI communication — the hardest, most error-prone part of a distributed FFT) that didn't need to be built from scratch.

**What is NOT reusable, and why**: `my_slab_based_fft_init()` only sets up layout bookkeeping (which task owns which slab) — it does **not** create the actual FFTW plan objects, despite superficially looking like a complete "init" function. Gravity's own plan-creation code (in `pm_nonperiodic.c`/`pm_periodic.c`) creates real-to-complex plans for the z-axis specifically (since gravity's density field starts real) and complex-to-complex for y/x (operating on the already-complex intermediate). Our wavefunction is complex *at every stage* — we need complex-to-complex plans for all three axes, which required writing our own plan-creation function, mirroring gravity's structure but genuinely different in content. **This exact gap (assuming `my_slab_based_fft_init` creates plans when it doesn't) caused a real segfault** the first time it was missed — see §3.6.

### 3.4 The Poisson solve — doubled mesh, and why gravity's own kernel doesn't apply

The isolated (non-periodic) potential update uses the standard Hockney & Eastwood zero-padding technique: a doubled mesh (`2×All.FDMGrid` per dimension), a Green's function kernel built once and cached, FFT-multiply-inverse-FFT each call. This deliberately reuses our own complex-to-complex machinery throughout (treating the real-valued density/kernel as complex with zero imaginary part) rather than adding a second, r2c-specific code path — costs 2× the memory/compute of a "proper" r2c implementation, but is dramatically simpler given everything already validated for c2c. A documented future optimization, not an oversight.

**Gravity's own non-periodic kernel (`pm_setup_nonperiodic_kernel` in `pm_nonperiodic.c`) cannot be reused directly, for two real reasons, not simplifications:**
1. It's a **long-range-only split kernel** (`fac = 1 - erfc(u)`, an Ewald/TreePM short-range/long-range split) — gravity's tree already handles short-range forces, so the PM mesh only needs the smoothed remainder. We have no separate tree handling the wavefunction's self-gravity; we need the **full, unsplit** `-1/r` potential.
2. It includes a **CIC deconvolution correction**, compensating for the Cloud-In-Cell interpolation error from depositing *particle* masses onto a grid. Our density (`|ψ|²`) is already exactly on the grid — applying this correction would actively introduce an error, not fix one.

Our own kernel uses the plain `-1/r` Green's function, with a simple half-cell softening at `r=0` (`-1/(0.5·dx)`) — explicitly flagged as an approximate, "basic version first" choice, not a precisely-derived analytic self-energy constant. Affects only a single cell's self-contribution; a secondary refinement if it ever matters.

### 3.5 Redistribution between two independently-decomposed FFT plans

`FDM_plan` (size `N`) and `FDM_PoissonPlan` (size `2N`) have **independently computed slab decompositions** — a task's local x-rows in one do not generally align with its local x-rows in the other, since they're different total sizes divided among the same `NTask`. There is no way to avoid explicit MPI communication between them by clever indexing alone. Both plans' `slab_to_task[]` arrays are fully populated on every task (via the `Allgather` inside `my_slab_based_fft_init`), so every task can independently compute its own send/receive plan without a separate negotiation round — this is what `fdm_redistribute_rows()` does.

### 3.6 Bugs found and fixed while building the Poisson solver — worth knowing as a class, not just individually

Five real bugs, in sequence, none caught by reasoning about the code alone — every one found by building an independent check (a from-scratch Python reference, a known analytic solution, or a direct diagnostic print) and comparing:

1. **Missing FFTW plan creation** for `FDM_PoissonPlan` — called `my_slab_based_fft_init` (layout only, see §3.3) and used the result directly without ever creating the actual plan objects. Segfault, caught immediately.
2. **LIFO free-order violation** — freed a scratch buffer mid-function while a later-allocated one was still alive; `mymalloc`'s own safety check caught it (same discipline learned the hard way during the SIDM work).
3. **Wrong physics scaling** — used `4πG` when just `G` was correct (the `-1/r` kernel already has the `4π` built into how it solves the Poisson equation), and was missing the `dx³` cell-volume factor the discrete FFT convolution needs to approximate the continuous integral it's standing in for. Found by building an independent Python implementation and comparing against the known analytic uniform-sphere potential — not by re-deriving the formula a second time and trusting it.
4. **Undersized MPI request buffer** — sized for "one request per task" when the actual pattern is one request per grid row; a genuine out-of-bounds write. Fixed, though it turned out not to be the cause of the remaining symptom (below) — worth knowing multiple real bugs can coexist and only some explain a given observed failure.
5. **The actual root cause of the largest remaining error**: source and destination row widths differ (`N` vs `2N`) between the two plans, and a flat contiguous row copy does not correctly place an `N×N` block into the right strided `(j,k)` positions of a `2N×2N` row. Needed an explicit pack/unpack rewrite (`fdm_redistribute_rows`'s current form), not a size fix.

**Validation after all fixes**: 1.3-2.2% relative error against the known analytic uniform-sphere potential, matching an independent Python reference to a fraction of a percent, identical results across 1/2/4 tasks.

### 3.7 Non-cosmological simplification (a=1)

The full paper's equations (17-21) include cosmological scale-factor (`a(t)`) dependence throughout. Since Phase 1's target is a non-cosmological central-galaxy test, all of `fdm_drift`, `fdm_kick`, `fdm_get_timestep`, and `fdm_update_potential` are implemented at `a=1` — a deliberate simplification for the current scope, not an oversight. Would need explicit revisiting before any cosmological (as opposed to isolated) FDM run.

### 3.8 `ℏ/m` and mass-scale unit conversions — cross-checked independently

`All.FDMMass` is `mc²` in eV (matching the literature's own convention for quoting FDM particle masses, not the mass itself). The `ℏ/m` code-units conversion (`fdm_hbar_over_m_code()`) was cross-checked against the paper's own quoted `λ_dB` value (1.21 kpc at `mc²=10⁻²² eV`, `v=100 km/s`) before being trusted — matched to 3 significant figures. The same "derive independently, cross-check against a citable external number" approach that validated SIDM's `σ/m` conversion.

**Worth flagging**: the actual numerical mass scale involved (`m_code ~ 10⁻⁹⁹` in these units, for an ultralight boson) is astronomically smaller than typical N-body particle masses. A synthetic test constructing `ψ` by inverting a desired density (`ψ = √(ρ/m)`) can produce extremely large intermediate `ψ` values (`~10⁴⁴` and up) — comfortably within double precision's range, but **will silently overflow to infinity in single precision** (max `~3.4×10³⁸`). This is exactly what happened when a test was accidentally run against the sandbox's default single-precision config rather than the project's actual `DOUBLEPRECISION_FFTW` build — see the practical notes below on always testing against the actual build configuration in use.

### 3.9 Phase 1 completion: wiring into `run.c`, and the soliton stability validation

`fdm_advance_to_time()` is called once per outer main-loop iteration, immediately after `find_next_sync_point()` updates `All.Ti_Current` — FDM runs on its own fully independent clock, sub-cycling internally via `fdm_step()` as many times as its own timestep criterion demands to land exactly on the target time, rather than participating in Arepo's hierarchical timebin system.

**Two real bugs found integrating this into an actual run, neither in FDM's own module code**:
1. **A movable-memory allocation-order bug**: `fdm_update_potential()`'s and `fdm_redistribute_rows()`'s scratch buffers were originally allocated fresh on every call (potentially hundreds of times per `fdm_advance_to_time()`). `mymalloc.c`'s own documentation states a movable block can only be safely shifted if all *subsequent* allocated blocks are also movable — Arepo's own particle arrays (`P[]`/`SphP[]`, movable, resized during domain decomposition) are allocated *after* FDM's arrays. Non-movable scratch buffers existing after those particle arrays whenever `fdm_update_potential()` happened to be running could corrupt the movable-block bookkeeping the moment a particle-array resize occurred at that exact instant. Fix: pre-allocate all scratch buffers once, early (`fdm_gradient_allocate()`-style pattern, established here first and reused throughout Phase 2a), before particle arrays exist at all.
2. **A genuine gap in Arepo's own base code**, not something this module introduced: `All.Asmth[0]`/`All.Rcut[0]` are only ever set in the *periodic* PM code path (`pm_periodic.c`) — `pm_nonperiodic.c` never sets them for the plain non-periodic case (only the `PLACEHIGHRESREGION` zoom-region values). Left at zero, this causes a division-by-zero propagating to an undefined integer cast and a crash the moment any particle pair's gravity gets evaluated, for *any* `PMGRID`+`GRAVITY_NOT_PERIODIC`+no-zoom-region run — not FDM-specific, just never previously exercised. Fixed by mirroring the periodic case's exact formula in `pm_nonperiodic.c`.

**Soliton IC generation and validation** (`make_soliton_ic.py`, not part of the C module): constructs a pure, isolated soliton from the Schive et al. (2014b) analytic profile (`ρ(r) = ρ₀/(1+0.091(r/rc)²)⁸`) — real-valued `ψ`, zero phase everywhere, matching the soliton's defining property as a stationary ground state. A **placeholder particle mechanism** (`make_placeholder_particles.py`) was needed alongside it: Arepo's framework assumes at least one particle exists per task (`test_id_uniqueness()` terminates on `NumPart==0`), and Phase 1 has no particles at all otherwise — minimal, negligible-mass Type-1 particles satisfy this without participating in any real physics.

**Multi-resolution validation (N=64/128/256)** found:
- **A genuine, consistent breathing-mode oscillation** in the core radius — same period and phase across all three resolutions (strong evidence of real physics, not a numerical artifact), with amplitude that shrinks with resolution at roughly the rate expected for genuine numerical convergence (a `~3.2×` reduction per doubling, close to the `4×` expected for second-order convergence), converging toward but not yet at a stable asymptotic value at N=256.
- **A small outer-envelope feature** (mass apparently shed into the outskirts over time) that looks essentially resolution-independent between N=128 and N=256 — arguing against pure discretization as the explanation, though this hasn't yet been disentangled from a possible box-size (boundary) artifact, since all runs used the same `L=20 kpc`.
- **Narrow, sharp dips in the outer-envelope density at specific, resolution-independent snapshot indices** — plausibly genuine wave-interference nodes (the wavefunction's amplitude passing near zero at specific times/locations), not investigated further but consistent with this project's own original interest in interference structure.

Analysis tooling (`fdm_field_tools.py`, `validate_soliton_stability.py`) was built and debugged against these real runs, not just written and assumed correct — see the transcript for the specific bugs caught this way (a units/mass-division error in the analytic overlay, comparing a coarse radial bin average against a steeply-peaked profile's true central value, and separately, a file-mixup investigation that turned out to be user error rather than a code bug, resolved by adding the actual grid size `N` directly into every plot title/filename as a structural safeguard against exactly that class of confusion recurring).

### 3.10 Phase 2a: why not `generic_comm_pattern` — a different, simpler communication pattern

Arepo's most commonly-reused inter-task communication machinery, `generic_comm_pattern` (`src/utils/generic_comm_helpers2.h`), is deeply tied to **tree-walk** infrastructure (`NodeList`, `Firstnode`, `NTopleaves`) — built for particles that need to *search* which tree nodes, local or remote, to visit for a neighbor-finding or force-summation problem. Mesh-to-star interpolation and star-to-mesh deposition have no search to do: a star's position directly and cheaply determines which task(s) own the relevant mesh slab, via `FDM_plan.slab_to_task[]`. Reusing `generic_comm_pattern` here would mean fabricating a node-list mechanism for a problem that doesn't have one.

Mirrored a different, simpler, already-proven pattern instead: gravity's own PM particle↔mesh communication (`pmforce_nonperiodic_uniform_optimized_prepare_density`/`readout_forces_or_potential` in `pm_nonperiodic.c`) — count which task(s) each particle's mesh footprint touches, exchange via a plain `myMPI_Alltoallv`, compute locally on whichever task owns each corner, and (for the mesh→star direction only) exchange results back using the same send/recv role swap. Since `FDM_plan` only distributes X (Y and Z are always fully local to every task), a star's trilinear stencil touches at most 2 tasks — simpler than gravity's own `FFT_COLUMN_BASED` case (up to 4 tasks, since that decomposition splits both X and Y); only gravity's simpler, non-column-based branch was needed as a reference.

**Deposition (star→mesh) is one-way**; interpolation (mesh→star) is a round trip. Verified this asymmetry directly from gravity's own code (`rhogrid[...] += mass*weight`, no corresponding "send back" step for density deposition) rather than assuming symmetry between the two directions.

### 3.11 The transpose layout bug — the most expensive lesson from this session, worth internalizing as a class

Computing `FDM_ForceX` (the force component requiring cross-task communication, since X may be distributed) needed a transpose so X becomes locally-owned. `fdm_transpose_x_y_A`/`B` (mirroring gravity's own `my_slab_transposeA`/`B`, but with stride `N` instead of `Ngrid2`, since this module's arrays are deliberately unpadded — see §3.3) were faithful, correct copies of gravity's proven algorithm. **The bug was in how the caller indexed into their output**, not in the transpose functions themselves: the post-transpose array has the *full-range* dimension outermost (largest stride) and the *local* dimension in the middle — confirmed by gravity's own `TI(x,y,z) = GRID*(x + y*nslab_x) + z` macro, the authoritative specification — and an earlier version of this code had these backwards.

**Why this went undetected initially, and the generalizable lesson**: the first validation used a *linear* test potential (`Φ = x - cx`), whose derivative is constant everywhere. A constant-slope test cannot distinguish correct indexing from swapped indexing — it gives the right answer either way, since a linear function's slope doesn't depend on *where* along the axis you sample it. It took switching to a genuinely position-dependent test (`Φ = x²+y²+z²`, whose gradient varies with position) to expose the bug — a relative error of `~2.0`, the exact signature of a sign flip, which led to fixing an *actual*, separate sign bug (see below) that didn't fix the real problem, before eventually tracing it to the layout swap via a long sequence of isolating tests (checking the transpose round-trip alone, a trivial same-position copy through the identical call structure, an intermediate-state check comparing single-task vs multi-task behavior — the single-task case being *independently* degenerate and misleading, since a single task's "transpose" reduces to a byte-for-byte identity copy that trivially "round-trips" regardless of whether the index *meaning* is right).

**A second, genuinely separate bug found in the same investigation**: `fdm_diff1d_strided`'s stencil already returns `-f'(x)` (i.e., the force itself, `-∇Φ`) directly — an earlier version applied an *additional* negation on top, exactly cancelling the correct sign. Caught by the same quadratic-potential test (relative error `~2.0`, before the layout bug was found underneath it) — a useful reminder that more than one real bug can coexist and manifest through the same symptom.

**Validated properly, eventually**: exact agreement with the analytic gradient to floating-point precision (`~10⁻¹⁵` absolute, `~10⁻¹⁷` relative) across 1, 2, 3, 4, and 8 tasks — including the 3-task (odd, asymmetric decomposition) case specifically, since powers of two alone would not have been a sufficiently different test of the general decomposition logic.

### 3.12 The star coupling itself: interpolation, deposition, and integration timing

**`fdm_interpolate_to_stars()`** (mesh→star): trilinear interpolation of `FDM_Potential` and all three force components together, in one communication round (not four separate ones — a star needs all four simultaneously, and the exchange overhead is the same either way).

**`fdm_deposit_star_mass()`** (star→mesh): CIC deposition of stellar mass into `FDM_StarMassDensity`, one-way (§3.10). This stores *deposited mass* per cell (CIC weights sum to exactly 1.0 per particle), not density — dividing by cell volume converts it to the same density units `FDM_RhoLocal` already uses, before adding it into the same source term the wavefunction's own density feeds (`fdm_update_potential()`, `fdm_poisson.c`).

**A real LIFO bug found in both functions** (same root cause, same fix, found once and then checked for and confirmed absent in the second occurrence): allocating `partout` before `partin`, then trying to free `partout` immediately after the exchange while `partin` (allocated after it) was still alive and needed. Gravity's own code allocates in the opposite order (`partin` first) specifically to avoid this — confirmed by re-reading gravity's own allocation order directly rather than assuming the "obvious" order was safe.

**Integration timing, and the limitation this represents**: `fdm_deposit_star_mass()` is called once per `fdm_advance_to_time()` invocation (i.e., once per outer sync point), before *any* potential update that call performs — including the bootstrap call, so even the very first solve includes the stellar contribution. `fdm_interpolate_to_stars()` is called once, at the end, after FDM's own sub-cycling has fully caught up to the target time. This is deliberate and correct *given* that star positions don't change during FDM's own finer internal sub-cycling (stars only move via Arepo's own gravity/timestep machinery, at the outer sync-point cadence) — but it is a genuine, explicitly-flagged **basic-version-first simplification**: stars feel FDM, and FDM feels stars, only at sync-point boundaries, not at any finer time resolution during the sub-cycling in between. Worth revisiting if the actual science ever needs finer coupling than that.

**Hooking the force into Arepo's own gravity accumulation** (`gravity_force_finalize()`, `accel.c`): added *after* `P[i].GravAccel[j] *= All.G`, not before. `FDM_ForceX/Y/Z` is `-∇Φ` where `Φ` is a potential per unit mass (this module's consistent convention throughout) — already an acceleration, already including `G` (baked into `Φ`'s own normalization in `fdm_poisson.c`) — adding it before the `*=All.G` step would double-apply `G` to this term specifically. Gated to `P[i].Type==4` (stars), matching Phase 2a's explicit scope.

**End-to-end validation, not just per-function unit tests**: a real simulation run (the soliton IC, 64 placeholder-turned-star particles, both `SIDM`-enabled and `SIDM`-fully-disabled configurations) completed 32 sync points cleanly with the full coupling active. Checked actual physics, not just absence of a crash: all 64 stars, which started at exactly zero velocity, picked up genuinely nonzero, star-to-star-varying velocities — direct evidence the mesh→star force is actually being felt, not silently no-op'ing. Confirmed byte-for-byte identical results with `SIDM` present vs. absent, closing out an open question about whether any accidental cross-module dependency had crept in (there wasn't one — zero references to `SIDM` anywhere in the FDM module or the `accel.c` hook).

---

## 4. File manifest

See `IMPLEMENTATION.md` for a complete, function-by-function reference. Summary:

```
src/fdm/
  fdm.h                     — module header: FDM_plan, FDM_Psi, FDM_PoissonPlan, FDM_Kernel, FDM_Potential,
                              FDM_ForceX/Y/Z, FDM_StarResult, FDM_StarMassDensity, all function prototypes,
                              PMGRID compile-time dependency check
  fdm_field.c               — field allocation, base FFT plan creation (genuinely complex-to-complex)
  fdm_integrator.c          — drift, kick, timestep criterion, the full kick-drift-kick step driver,
                              fdm_advance_to_time() (the run.c-facing entry point, now also driving the
                              Phase 2a star deposit/interpolate calls)
  fdm_poisson.c              — doubled-mesh Poisson solve: kernel construction, potential update (now
                              including the stellar source term), cross-plan row redistribution
  fdm_gradient.c             — Phase 2a: force-on-mesh via 4th-order finite differencing, including the
                              transpose-based X-direction handling (§3.11)
  fdm_particle_coupling.c    — Phase 2a: fdm_interpolate_to_stars(), fdm_deposit_star_mass(), FDM_StarResult
                              and FDM_StarMassDensity's own allocation

  fdm_test_fft.c             — round-trip + Parseval validation of the base FFT machinery
  fdm_test_single_mode.c     — validates the exact post-transform index-to-k-vector mapping
                              (needed before trusting the drift operator's k^2 computation)
  fdm_test_drift.c           — norm conservation + exact analytic phase prediction for the drift operator
  fdm_test_kick.c            — pointwise magnitude preservation + exact phase prediction for the kick operator
  fdm_test_poisson.c         — validates the Poisson solve against a known analytic uniform-sphere potential
  fdm_test_step.c            — validates the full kick-drift-kick sequence (norm conservation over multiple steps)
  fdm_test_force.c           — validates fdm_compute_force() against a known quadratic potential (§3.11)
  fdm_test_interpolation.c   — validates fdm_interpolate_to_stars() (analytic trilinear-interpolation error bound)
  fdm_test_deposit.c         — validates fdm_deposit_star_mass() (mass conservation + grid-aligned exact case)

(not in src/fdm/, but part of this module's own tooling)
  make_soliton_ic.py              — constructs the Schive et al. soliton IC used for validation
  make_placeholder_particles.py   — minimal Type-1 placeholder particles (Phase 1, before any real coupling)
  make_star_particles.py          — minimal Type-4 star placeholders (Phase 2a end-to-end testing)
  fdm_field_tools.py               — FDMFieldTools class for reading FDM field snapshots (mirrors SnapshotTools)
  validate_soliton_stability.py    — the soliton stability validation script (density profile, core radius,
                                     norm conservation over time, resolution comparisons)
```

Test files are standalone executables (their own `main()`, linked against the main build's object files with `main.o` swapped out) — not part of the production Arepo binary, and not yet wired into any automated test harness.

---

## 5. Required reading

1. **This file**, then `IMPLEMENTATION.md` for the precise code reference.
2. `src/fdm/fdm.h` — the PMGRID dependency and why (§3.3), and every array/plan's role.
3. `src/fdm/fdm_poisson.c`'s header comments — the doubled-mesh technique, why gravity's kernel doesn't apply (§3.4), and the redistribution mechanism (§3.5) — the most involved and most-corrected piece of this module.
4. `src/fdm/fdm_gradient.c`'s header comments — the transpose layout derivation (§3.11), the single most hard-won piece of Phase 2a.
5. May & Springel (2021), "Structure formation in large-volume cosmological simulations of fuzzy dark matter" — Eqs. 4-23 specifically for the theoretical background and the split-step method this module implements (with the reuse/departures documented in §3 above).
6. Schive et al. (2014b) — the soliton analytic profile and core-halo relation used in §3.9's validation.
7. The BECDM/"Galaxy Formation with wave/fuzzy dark matter" paper series and "Fuzzy Gasoline" (2024) — for context on what's already been done elsewhere (§1), before extending this module's science scope further.

---

## 6. TODO / known limitations, roughly in priority order

### Resolved since Phase 1's original writing
- [x] ~~Not yet wired into Arepo's actual run loop~~ — done, see §3.9.
- [x] ~~No initial-condition generation for a real target halo~~ — the pure soliton case is done (§3.9); the fuller soliton+NFW-envelope halo construction (needing the velocity-field phase-solving machinery, Eqs. 26-27) remains genuinely deferred, see below.
- [x] ~~Baryonic coupling (Phase 2a)~~ — done, see §3.10-3.12.

### Blocking "Phase 2a is actually usable for real science", not just structurally complete
- [ ] **Sync-point-only coupling** (§3.12) — stars feel FDM, and vice versa, only at the outer sync-point cadence, not during FDM's own finer sub-cycling. Explicitly flagged as basic-version-first; would need revisiting if the science requires finer time resolution than that.
- [ ] **Y/Z boundary handling for stars near the box edge is clamped, not properly handled** — `fdm_trilinear_accumulate()` and the deposit function's own corner logic both clamp `yy`/`zz` at the box edge rather than doing anything more careful. Fine given this project's stars are expected well within the box (matching the existing central-galaxy scoping), but worth revisiting if that assumption changes.
- [ ] **The pure-soliton IC's outer-envelope feature and breathing-mode amplitude are not yet fully converged** (§3.9) — worth another resolution step or a dedicated box-size test (to disentangle genuine mass-shedding from a boundary artifact) before treating the soliton validation as fully closed out, independent of Phase 2a itself.

### Architecture, deferred by deliberate choice
- [ ] Gas coupling (Phase 2b) — not started; the genuinely hard part, per §2.
- [ ] BH coupling (Phase 2c) — not started.
- [ ] The fuller soliton+NFW-envelope halo IC (needs the velocity-field phase-solving machinery, Eqs. 26-27) — deferred; the pure soliton was deliberately scoped as the first, simpler validation target instead (§3.9).
- [ ] Hybrid zoom capability (Phase 5) — scoped, not started; see §2. Phase 2a's transpose-based, task-ownership-driven communication pattern (§3.10) is directly relevant groundwork here.
- [ ] `r2c` optimization for the Poisson solve — currently uses c2c throughout (2× memory/compute cost), deliberately deferred (§3.4).
- [ ] Cosmological (`a≠1`) support — not implemented (§3.7).
- [ ] `FFT_COLUMN_BASED` — not handled; `fdm_field.c`/`fdm_poisson.c` both `#error` out if this is set. Only the slab-based path is implemented (and the Phase 2a communication pattern, §3.10, is only implemented for gravity's simpler, non-column-based branch). Would matter for very large task counts (the column-based decomposition scales better, per `pm_mpi_fft.c`'s own two-strategy design), not for current-scale testing.
- [ ] Consecutive-half-kick merging optimization (a performance optimization the paper itself describes, combining the trailing half-kick of one step with the leading half-kick of the next) — deliberately not implemented, correctness-neutral, deferred.

### Smaller, worth tracking
- [ ] `r=0` kernel self-term uses a simple half-cell softening, not a precisely-derived analytic constant (§3.4).
- [ ] No automated test harness — all test executables are run manually, against whatever `Config.sh` happens to be active; easy to accidentally test against the wrong precision/flag configuration (see practical notes below).
- [ ] FDM output cadence is still every outer sync point (a "basic version first" simplification flagged since Phase 1), not integrated with Arepo's own output-time-list mechanism — may write far more often than actually wanted for a real run.

---

## 7. Practical notes for whoever runs this next

- **Always verify which `Config.sh` is actually active before trusting a test result**, especially `DOUBLEPRECISION_FFTW`. A test run against the sandbox's default (single-precision) config produced a wall of `NaN` that took real effort to diagnose, before realizing it was single-precision overflow in a value that's completely benign in double precision — not a logic bug at all. This project's actual target build (`DOUBLEPRECISION_FFTW` on) should be the one every test is ultimately checked against, even if faster iteration happens against a different config first.
- **When two independently-sized/decomposed data structures need to exchange data (as with the two FFT plans here), do not assume a flat/contiguous copy is correct** — verify the actual row/stride structure on both sides explicitly. This bug (§3.6, item 5) was the single largest source of wasted effort in Phase 1.
- **Build independent verification before trusting a new physics formula** — the `4πG`/`dx³` scaling bug (§3.6, item 3) was only caught by writing a from-scratch Python implementation and comparing against a known analytic solution, not by re-checking the derivation on paper a second time. This is the same lesson that validated SIDM's cross-section units and this module's own `ℏ/m` conversion — it generalizes.
- **A test potential with a constant derivative (linear in position) cannot detect position-dependent indexing bugs** — it will pass even when the indexing is genuinely wrong, since the "right answer" doesn't depend on where you sample it (§3.11). Always validate against a function whose derivative actually varies with position (a quadratic is the simplest choice) before trusting a gradient/force computation, even if a linear test already "passed."
- **When reusing gravity's own communication patterns, check the authoritative macro/formula directly rather than infer the layout from testing alone** — an empirical test on a single task can be *independently* misleading here, since a single task's "transpose" degenerates to a byte-for-byte identity copy, which trivially "passes" a round-trip test regardless of whether the index *meaning* on either side is actually correct (§3.11). Gravity's own `TI`/`FI`-style macros, where they exist, are the ground truth — read them directly rather than reverse-engineer the convention from behavior.
- **When allocating a matched send/receive buffer pair that you intend to free one of immediately after the exchange, check which one gravity's own code allocates first** — the convention (allocate the one you'll free immediately *last*) exists specifically so that free is LIFO-safe; getting this backwards produces an immediate, easy-to-catch `mymalloc` safety-check failure, but only if you actually run the code rather than assume the "obvious" order is fine (§3.12).
- **A standalone test that creates identical particle data must still account for which task "owns" each particle**, if it's meant to mimic a real domain-decomposed run — every task independently creating the same test particles will silently multiply results by `NTask`, a bug in the test, not the code being tested, but one that produces a very clean, misleadingly systematic-looking failure pattern.
- Every multi-task bug in this module so far was, in fact, also visible on a single task *except* one: the Phase 2a transpose-layout bug (§3.11) specifically required a 2-task (or higher) test to disambiguate from a misleading single-task degeneracy. Confirming identical behaviour across a genuinely varied set of task counts (including an odd number, not just powers of two) remains valuable evidence that a fix is correct for the general case, not coincidentally correct for one specific decomposition — this session's 3-task test of the force computation is a concrete example of that mattering, not just a formality.
