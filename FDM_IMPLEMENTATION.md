# FDM Implementation Reference

Companion to `FDM.md`. That file explains architecture decisions and *why*; this file is a precise map of *what exists and where* — every struct, function, and global, so anyone can navigate the codebase without re-reading narrative. Not a substitute for reading the actual header/function comments in the source, which go into more depth on the trickier pieces (especially in `fdm_poisson.c`).

---

## `src/fdm/fdm.h`

Module header. Included by every FDM-touching file.

**Compile-time dependency check**: `#ifndef PMGRID / #error` — requires `PMGRID` set in `Config.sh`. Not because FDM uses gravity's own long-range PM solve; `fft_real`/`fft_complex`/`fft_plan` (in `dtypes.h`) and the FFTW3 library linkage are gated behind `PMGRID`, and this module reuses those type definitions rather than duplicating them.

**Globals — Phase 1, field/Poisson solve**:
- `fft_plan FDM_plan` — the base mesh (`All.FDMGrid` per dimension), used for `ψ`'s own drift/kick.
- `fft_complex *FDM_Psi` — the wavefunction field itself, distributed per `FDM_plan`'s slab decomposition. No real-to-complex padding anywhere — genuinely complex throughout, never a real-valued intermediate stage.
- `fft_complex *FDM_PsiWorkspace` — scratch space for `FDM_plan`'s FFT executions.
- `fft_real *FDM_Potential` — the gravitational potential, real-valued, same real-space distributed layout as `FDM_Psi`. Consumed by `fdm_kick`; produced by `fdm_update_potential`.
- `fft_plan FDM_PoissonPlan` — a *separate*, doubled-mesh (`2×All.FDMGrid` per dimension) plan, used only for the Poisson solve (Hockney & Eastwood zero-padding). Independently decomposed from `FDM_plan` — see `fdm_poisson.c`'s redistribution logic.
- `fft_complex *FDM_Kernel` — the FFT'd Green's function kernel on the doubled mesh, built once (`fdm_poisson_kernel_init`), cached, reused every potential-update call.
- `fft_complex *FDM_PoissonWorkspace` — scratch for `FDM_PoissonPlan`'s FFT executions.
- `fft_complex *FDM_RhoLocal`, `FDM_RhoPadded`, `FDM_PhiLocal`, `FDM_RedistPackBuf`, `void *FDM_RedistRequests` — persistent scratch buffers for `fdm_update_potential()`/`fdm_redistribute_rows()`, pre-allocated once in `fdm_poisson_kernel_init()` (called from `fdm_allocate()`, before particle arrays exist) rather than per-call. Not a performance choice — fixes a genuine movable-memory-block crash (`FDM.md` §3.9, item 1); this pattern (allocate everything mesh-sized once, early, before `P[]`/`SphP[]` exist) is reused throughout the rest of this module.

**Globals — Phase 2a, force computation** (`fdm_gradient.c`):
- `fft_real *FDM_ForceX`, `FDM_ForceY`, `FDM_ForceZ` — `-∇Φ`, same real-space layout as `FDM_Potential`. Computed by `fdm_compute_force()`.

**Globals — Phase 2a, particle coupling** (`fdm_particle_coupling.c`):
- `fdm_star_result` (typedef'd struct: `Potential`, `ForceX`, `ForceY`, `ForceZ`) and `fdm_star_result *FDM_StarResult` — per-*particle* (not mesh) results from `fdm_interpolate_to_stars()`. Genuinely different allocation pattern from everything else in this module: sized to `All.MaxPart`, grown in lockstep with `P[]` itself via `reallocate_memory_maxpart()` (`src/utils/allocate.c`), mirroring SIDM's `DMSP[]` array exactly rather than the mesh-array pattern above.
- `fft_real *FDM_StarMassDensity` — stellar mass deposited by `fdm_deposit_star_mass()`, same real-space distributed layout as `FDM_Potential`. Read directly by `fdm_update_potential()` (`fdm_poisson.c`) to add the stellar source term — a genuine, necessary cross-file dependency (not `static`, unlike some of this module's other internal scratch buffers).

**Prototypes**: `fdm_allocate`, `fdm_free`, `fdm_hbar_over_m_code`, `fdm_drift`, `fdm_kick`, `fdm_get_timestep`, `fdm_step`, `fdm_advance_to_time`, `fdm_poisson_kernel_init`, `fdm_update_potential`, `fdm_gradient_allocate`, `fdm_gradient_free`, `fdm_compute_force`, `fdm_particle_coupling_allocate`, `fdm_particle_coupling_free`, `fdm_interpolate_to_stars`, `fdm_deposit_star_mass`.

**New `All.` parameters** (added to `allvars.h`, not this file): `All.FDMGrid` (int, mesh size `N`), `All.FDMBoxSize` (double, physical box size `L`, code length units), `All.FDMMass` (double, `mc²` in eV — the boson rest-mass energy, not the mass itself), `All.FDMICFile` (string, path to the soliton/field IC file).

---

## `src/fdm/fdm_field.c`

Field allocation and the base FFT plan's actual creation.

**`void fdm_allocate(void)`** — the driver:
1. `my_slab_based_fft_init(&FDM_plan, N, N, N)` — layout bookkeeping only (slab assignment across tasks); does **not** create FFTW plan objects (a genuine gap in what this function does, easy to assume otherwise — see `FDM.md` §3.3).
2. Allocates `FDM_Psi`, `FDM_PsiWorkspace`, `FDM_Potential` (in that order — matters for the LIFO free discipline in `fdm_free`).
3. Creates the actual FFTW plan objects for all six transforms (`forward`/`backward` × `x`/`y`/`z` direction), using `FFTW(plan_many_dft)` (genuinely complex-to-complex) for **all three axes**, including `z` — this is the one deliberate departure from gravity's own template in `pm_nonperiodic.c`, which uses `plan_many_dft_r2c`/`c2r` for `z` specifically since gravity's density starts real.
4. Calls `fdm_poisson_kernel_init()` — which itself now also allocates `FDM_RhoLocal`/`FDM_RhoPadded`/`FDM_PhiLocal`/`FDM_RedistPackBuf`/`FDM_RedistRequests` (`fdm_poisson.c`).
5. Calls `fdm_gradient_allocate()` (`fdm_gradient.c`) — `FDM_ForceX/Y/Z` and the transpose scratch buffers.
6. Calls `fdm_particle_coupling_allocate()` (`fdm_particle_coupling.c`) — `FDM_StarMassDensity` (mesh-sized; `FDM_StarResult`, being per-particle, is allocated separately, in `src/utils/allocate.c`, alongside `P[]`/`SphP[]`).

Each of these three calls happens last-in-first-freed order relative to the others — see `fdm_free()` below.

**`void fdm_free(void)`** — frees everything in exact reverse allocation order: `fdm_particle_coupling_free()` first (most recently allocated), then `fdm_gradient_free()`, then the Poisson-related allocations (`FDM_RedistRequests` through `FDM_PoissonPlan`'s five bookkeeping arrays), then the base allocations (`FDM_Potential`, `FDM_PsiWorkspace`, `FDM_Psi`, `FDM_plan`'s five bookkeeping arrays).

---

## `src/fdm/fdm_integrator.c`

Drift, kick, timestep criterion, and the full step driver — everything except the Poisson solve.

**`double fdm_hbar_over_m_code(void)`** — converts `All.FDMMass` (eV, `mc²`) to `ℏ/m` in code units (`code_length²/code_time`), via CGS. Cross-checked against the source paper's own quoted `λ_dB` value before being trusted (see `FDM.md` §3.8). `PLANCK` (allvars.h) is `h`, not `ℏ=h/(2π)` — divided out here.

**`static inline double fdm_k_of_index(int idx, int N, double L)`** — converts a grid index (`0..N-1`) to its physical wavenumber, standard FFT frequency convention (negative frequencies wrap around above the Nyquist index).

**`void fdm_drift(double dt)`** — the drift step (Eq. 19's central term), entirely in Fourier space: forward FFT, multiply each mode by `exp(-i·(ℏ/m)·dt·k²/2)`, inverse FFT with `1/N³` normalization.

**Index-to-`(kx,ky,kz)` mapping used here** — `index = i_local_y·N² + kx·N + kz`, `ky_global = FDM_plan.slabstart_y + i_local_y`. This was derived by tracing `my_slab_transpose` by hand, then **independently verified** with a single-Fourier-mode test (`fdm_test_single_mode.c`) *before* being trusted in this function — not assumed correct from the derivation alone.

**`void fdm_kick(double dt)`** — the kick step: pointwise phase rotation of `ψ` by the local potential, `exp(-i·(m/ℏ)·dt·Φ(x))`. Real-space only, no FFT — much lower risk than drift, with a correspondingly stronger available invariant (pointwise `|ψ|` preservation, not just aggregate norm conservation).

**`double fdm_get_timestep(void)`** — Eq. 21's criterion (non-cosmological, `a=1`): `dt < min[(4/3π)·(ℏ/m)·Δx², 2π·(ℏ/m)/|Φ_max|]`. `|Φ_max|` is MPI-reduced globally across the whole distributed `FDM_Potential` array. Requires `FDM_Potential` already populated; returns just the drift term if no potential has been computed yet (`Φ_max=0`).

**`void fdm_step(double dt)`** — the full sequence: `fdm_kick(dt/2)` → `fdm_drift(dt)` → `fdm_update_potential()` → `fdm_kick(dt/2)`. Deliberately does **not** implement the paper's consecutive-half-kick merging optimization across repeated calls — a documented, correctness-neutral, deferred performance optimization. **Precondition**: `FDM_Potential` must already be valid before the *first* call in a simulation (the opening half-kick uses whatever potential currently exists) — every subsequent call is self-sufficient, since the potential gets refreshed before the closing half-kick of each step.

**`void fdm_advance_to_time(integertime target_ti)`** — the `run.c`-facing entry point, called once per outer main-loop iteration, immediately after `find_next_sync_point()`. FDM runs on its own fully independent clock (`FDM_CurrentTime`), sub-cycling via `fdm_step()` as many times as needed to land exactly on the target time (never overshooting), rather than participating in Arepo's hierarchical timebin system.
1. **Phase 2a**: calls `fdm_deposit_star_mass()` once, unconditionally, before anything else in this function — including before the bootstrap `fdm_update_potential()` call on the very first invocation, so even the first potential solve already includes the stellar source term.
2. Bootstraps on first call: sets `FDM_CurrentTime = All.TimeBegin`, computes the initial potential.
3. Sub-cycles `fdm_step()` until `FDM_CurrentTime` reaches `target_time`.
4. **Phase 2a**: calls `fdm_interpolate_to_stars()` once, after sub-cycling completes — populating `FDM_StarResult` with the now-current potential/force, ready for `accel.c`'s `gravity_force_finalize()` to pick up before the next kick.
5. Writes an FDM field snapshot (`fdm_write_field`), tagged by a monotonic counter — every call, not integrated with Arepo's own output-time-list mechanism (a known, flagged simplification, `FDM.md` §6).

The single deposit/interpolate call per invocation (not per inner `fdm_step()`) is deliberate and correct *given* that star positions don't change during this function's own sub-cycling (§3.12 of `FDM.md`) — but is the specific mechanism behind the "sync-point-only coupling" limitation documented there.

---

## `src/fdm/fdm_poisson.c`

The isolated (non-periodic) Poisson solve. The most involved, most-corrected file in this module — see `FDM.md` §3.4-3.6 for the full reasoning and bug history; this section is the terse reference only.

**`static void fdm_redistribute_rows(fft_complex *src_rows, fft_plan *src_plan, int src_row_width, fft_complex *dst_rows, fft_plan *dst_plan, int dst_row_width, int transfer_width)`** — moves the `[0,N)×[0,N)` `(y,z)` sub-block of each row between two independently-decomposed plans. Explicit parameters for source/destination row width, since `FDM_plan`'s rows (`N` wide) and `FDM_PoissonPlan`'s rows (`2N` wide) are genuinely different — a flat contiguous copy is **only** valid when both sides have equal row width, which is not generally true here. Packs into a contiguous scratch buffer before sending when the source row is wider than `transfer_width`; unpacks a contiguous received buffer into the correctly-strided destination positions when the destination row is wider. Uses `gx` (global x-index) as the MPI tag throughout, so sends/receives correctly match regardless of which plan is playing the "source" vs "destination" role in a given call (the function is used symmetrically, once forward and once backward, per potential-update call).

**`void fdm_poisson_kernel_init(void)`** — one-time setup, called from `fdm_allocate()`:
1. `my_slab_based_fft_init(&FDM_PoissonPlan, 2N, 2N, 2N)`.
2. Allocates `FDM_Kernel`, `FDM_PoissonWorkspace`.
3. **Creates the actual FFTW plan objects** for `FDM_PoissonPlan` (this step was missing entirely in an earlier version — see `FDM.md` §3.6, item 1 — caused a segfault, fixed).
4. Builds the real-space kernel: plain `-1/r` Green's function (not gravity's erfc-split, CIC-deconvolved kernel — see `FDM.md` §3.4 for why), using the same wraparound indexing convention as gravity's own kernel construction. `r=0` self-term: simple half-cell softening (`-1/(0.5·dx)`), explicitly flagged as approximate.
5. FFTs the kernel once; result overwrites `FDM_Kernel` in place.

**`void fdm_update_potential(void)`** — the per-call potential update (Eq. 20c), no mean-density subtraction (that correction is specific to the periodic/cosmological case; this solves the true non-periodic equation directly, matching gravity's own `GRAVITY_NOT_PERIODIC` choice):
1. Computes `ρ = m·|ψ|²` (Eq. 5) on `FDM_plan`'s local rows — **no `4πG` factor here** (a bug in an earlier version — the `-1/r` kernel already has an implicit `4π` normalization built in; including it again double-counts the `4π` — see `FDM.md` §3.6, item 3). **Phase 2a**: adds `FDM_StarMassDensity[idx] / cell_volume` to this same source term — `FDM_StarMassDensity` stores CIC-*deposited mass* per cell (weights sum to exactly 1.0 per particle), not density; dividing by cell volume converts it to the same density units this term already uses.
2. Redistributes into the padded mesh (`fdm_redistribute_rows`, forward direction) — the destination buffer is zeroed first, which is what actually implements the zero-padding.
3. FFT, multiply by `FDM_Kernel`, inverse FFT, `1/(2N)³` normalization.
4. Redistributes back (`fdm_redistribute_rows`, reverse direction — the same function, src/dst roles swapped).
5. Extracts into `FDM_Potential`, applying `G_Newton · dx³` — the correct scaling, derived and verified numerically against the known analytic uniform-sphere solution (not from memory) — see `FDM.md` §3.6, item 3, for the missing-`dx³` half of this same bug.

---

## `src/fdm/fdm_gradient.c` (Phase 2a)

Computes `FDM_ForceX/Y/Z = -∇(FDM_Potential)` via 4th-order finite differencing. See `FDM.md` §3.11 for the full, hard-won derivation and bug history — this section is the terse reference only.

**`static inline double fdm_diff1d_strided(fft_real *field, size_t base, int idx, int N, size_t stride, double dx)`** — the verified stencil: `(4/3)(f[i-1]-f[i+1]) - (1/6)(f[i-2]-f[i+2])`, divided by `2·dx`, for interior points (`2 <= idx <= N-3`); falls back to a 2-point centered difference one cell from the boundary, and a 1-point one-sided difference at the absolute edge (`idx==0` or `N-1`) — deliberately not wrapped around the array boundary, since this is an isolated (non-periodic) domain, not a periodic one. Verified via symbolic Taylor expansion (not assumed) that the raw stencil equals `-2·dx·f'(x) + O(dx⁵)` — the function already returns the *force* (`-f'(x)`) directly, not the raw derivative; no additional sign flip needed at the call site (a real bug in an earlier version did exactly that).

**`static void fdm_transpose_x_y_A/B(fft_plan *plan, fft_real *field, fft_real *scratch)`** — faithful copies of gravity's own `my_slab_transposeA`/`B` (`pm_mpi_fft.c`) algorithm and MPI communication pattern, with stride fixed to `N` instead of `Ngrid2` (gravity's r2c-padding-specific stride, which doesn't apply here — this module is genuinely unpadded throughout, `FDM.md` §3.3). **Post-transpose layout** (the critical, hard-won detail, `FDM.md` §3.11): the full-range (original X) dimension is outermost (largest stride), the local Y-slab dimension is the middle dimension — confirmed via gravity's own `TI(x,y,z)` macro, not inferred from testing.

**`void fdm_compute_force(void)`**:
1. Z (stride 1) and Y (stride `N`) differencing: both purely local, no transpose needed (a task holds all Y/Z data for its own X-slabs).
2. X differencing: copies `FDM_Potential` into `FDM_TransposedPotential`, transposes via `fdm_transpose_x_y_A`, differences along the now-local former-X dimension (same stencil, same code, different stride/position in the loop nesting per the layout above), transposes the result back via `fdm_transpose_x_y_B`, copies into `FDM_ForceX`.

**`void fdm_gradient_allocate/free(void)`** — allocates/frees `FDM_ForceX/Y/Z` and this file's three dedicated scratch buffers (`FDM_TransposedPotential`, `FDM_TransposedForceX`, `FDM_CommScratch` — the latter used internally by the transpose functions' own communication). Called from `fdm_allocate()`/`fdm_free()`. `FDM_TransposedPotential`/`FDM_TransposedForceX` are sized using `FDM_plan.largest_y_slab` (the *global* maximum y-slab size across all tasks) — every task must size its own buffer for whatever it might receive after transposing, not just its own share.

---

## `src/fdm/fdm_particle_coupling.c` (Phase 2a)

The mesh↔star coupling: `fdm_interpolate_to_stars()` (mesh→star, a round trip) and `fdm_deposit_star_mass()` (star→mesh, one-way). See `FDM.md` §3.10-3.12 for the full architecture reasoning (why not `generic_comm_pattern`, the gravity-PM-mirroring communication pattern, the LIFO bug, integration timing).

**Communication pattern, common to both functions**: (1) count which task(s) each star's mesh footprint touches, via `FDM_plan.slab_to_task[]` on `slab_x`/`slab_x+1` (at most 2 tasks, since only X is distributed); (2) exchange positions (interpolation) or positions+mass (deposition) via a plain `myMPI_Alltoallv`; (3) compute/deposit locally on whichever task owns each corner; (4, interpolation only) exchange results back using the same send/recv role swap, reassembled using the *identical* iteration order as steps 1-2.

**`static void fdm_trilinear_accumulate(...)`** — the 8-corner trilinear weighting, mirroring gravity's own readout formula (`pm_nonperiodic.c`) exactly, accumulating `Potential`/`ForceX/Y/Z` together from a single task's local mesh data. Y/Z boundary handling is clamped (not wrapped), not fully general — see `FDM.md` §6 for the known limitation this represents.

**`static inline int fdm_particle_in_scope(int i)`** — `P[i].Type == 4` (stars). The single point of control for Phase 2a's particle-type scope; would need extending here first if gas/BH coupling (2b/2c) ever needs the same communication pattern for a different particle type.

**`void fdm_interpolate_to_stars(void)`** — populates `FDM_StarResult[i]` for every in-scope, in-box particle. Uses `struct fdm_star_partbuf` (position only) for the outbound exchange and `fdm_star_result` (4 doubles: potential + 3 force components) for the return — one combined round trip, not four separate ones.

**`void fdm_deposit_star_mass(void)`** — zeroes `FDM_StarMassDensity` at the start of every call (an accumulator, not a running total across calls), then deposits via `struct fdm_star_massbuf` (position + mass). One-way: no results exchanged back, matching gravity's own `rhogrid[...] += mass*weight` deposition exactly (verified directly from gravity's code, not assumed).

**`void fdm_particle_coupling_allocate/free(void)`** — allocates/frees `FDM_StarMassDensity` only (mesh-sized, so it follows this module's usual pre-allocate-early pattern). `FDM_StarResult` is per-*particle*, not mesh data — it's allocated separately, in `src/utils/allocate.c`, mirroring `DMSP[]`'s exact pattern (sized to `All.MaxPart`, grown via `reallocate_memory_maxpart()`, no independent reallocate function of its own needed since it's always tied 1:1 to `All.MaxPart`).

**Test-only debug helpers** (`fdm_debug_sum_star_mass_density`, `fdm_debug_get_star_mass_density_at`, `fdm_debug_get_transposed_force_x` in `fdm_gradient.c`, `fdm_debug_transpose_roundtrip` in `fdm_gradient.c`) — expose otherwise-private module state for validation; not part of the production code path, left in place rather than stripped out, since they're cheap and directly useful if this module's communication/layout logic ever needs re-diagnosing.

---

## `src/gravity/accel.c` (modified, not new — Phase 2a's hook into Arepo's own gravity)

**`gravity_force_finalize(int timebin)`** — modified to add `FDM_StarResult[i].ForceX/Y/Z` into `P[i].GravAccel[j]`, gated to `P[i].Type==4`, immediately *after* the existing `P[i].GravAccel[j] *= All.G` loop, not before. `FDM_ForceX/Y/Z` is `-∇Φ` where `Φ` is a potential per unit mass (this module's consistent convention) — already an acceleration, already including `G` (baked into `Φ`'s own normalization in `fdm_poisson.c`) — adding it before the `*=All.G` step would double-apply `G` to this term specifically. The only file outside `src/fdm/` this module modifies for Phase 2a (besides `src/utils/allocate.c`'s `FDM_StarResult` allocation, mirroring `DMSP[]`).

---

## Test files — standalone executables, not part of the production binary

Each has its own `main()`; built by compiling the test `.c` file, then linking against the main build's object files with `main.o` excluded. Not wired into any automated harness — run manually, and **must be checked against the actual `Config.sh` in use** (see `FDM.md` §7 on the single/double-precision FFTW pitfall).

**These are Phase 2a's actual validation suite** — the direct answer to "what validates this phase": `fdm_test_force.c`, `fdm_test_interpolation.c`, `fdm_test_deposit.c` below, plus the end-to-end run described in `FDM.md` §3.12. All three unit tests deliberately use a **quadratic**, not linear, test potential where position-dependence matters (`FDM.md` §3.11) — a linear test's constant derivative cannot detect a position-dependent indexing bug, since the "right answer" doesn't depend on where you sample it.

| File | Validates | Method |
|---|---|---|
| `fdm_test_fft.c` | Base FFT machinery | Round-trip recovery + Parseval's theorem (both index-order-agnostic) |
| `fdm_test_single_mode.c` | Post-transform index→`(kx,ky,kz)` mapping | Single known Fourier mode, checks the spike lands exactly where derived |
| `fdm_test_drift.c` | Drift operator | Norm conservation + exact analytic phase for a known mode |
| `fdm_test_kick.c` | Kick operator | Pointwise `|ψ|` preservation + exact phase, spatially-varying potential |
| `fdm_test_poisson.c` | Poisson solve | Known analytic uniform-sphere potential (interior + exterior) |
| `fdm_test_step.c` | Full `fdm_step` sequence | Norm conservation over multiple steps |
| `fdm_test_force.c` (Phase 2a) | `fdm_compute_force()` | Quadratic potential (`Φ=x²+y²+z²`), exact analytic gradient — floating-point-precision agreement in the interior, checked on 1/2/3/4/8 tasks (the 3-task case specifically needed to catch the transpose-layout bug, `FDM.md` §3.11) |
| `fdm_test_interpolation.c` (Phase 2a) | `fdm_interpolate_to_stars()` | Same quadratic potential; potential error checked against the *analytic* trilinear-interpolation error bound (not floating-point precision — trilinear interpolation is genuinely not exact for a quadratic); force error checked to be exact (force is linear in position, and trilinear interpolation *is* exact for linear functions) |
| `fdm_test_deposit.c` (Phase 2a) | `fdm_deposit_star_mass()` | Two checks: total mass conservation (CIC weights sum to exactly 1.0 per particle — floating-point-precision agreement, a strong position-independent property) and a grid-aligned star depositing its entire mass into exactly one cell, zero elsewhere (an exact, not approximate, check) |

---

## Quick lookup: "where is X actually computed?"

| Quantity | Computed in | Stored in |
|---|---|---|
| `ψ` evolution (kinetic term) | `fdm_drift()` | `FDM_Psi` (in place) |
| `ψ` evolution (potential term) | `fdm_kick()` | `FDM_Psi` (in place) |
| `Φ` (gravitational potential) | `fdm_update_potential()` | `FDM_Potential` |
| Green's function kernel | `fdm_poisson_kernel_init()` (once) | `FDM_Kernel` |
| Timestep | `fdm_get_timestep()` | returned, not stored |
| `ℏ/m` (code units) | `fdm_hbar_over_m_code()` | returned, not stored (cheap, recomputed each call) |
| Row redistribution between the two plans | `fdm_redistribute_rows()` | writes into whichever `dst_rows` buffer the caller passes |
| Force field (`-∇Φ`) (Phase 2a) | `fdm_compute_force()` | `FDM_ForceX`/`FDM_ForceY`/`FDM_ForceZ` |
| Mesh→star interpolated potential/force (Phase 2a) | `fdm_interpolate_to_stars()` | `FDM_StarResult[i]` (per-particle) |
| Star→mesh deposited mass (Phase 2a) | `fdm_deposit_star_mass()` | `FDM_StarMassDensity` (added into the source term by `fdm_update_potential()`) |
| Star's total gravitational acceleration, including FDM (Phase 2a) | `gravity_force_finalize()` (`accel.c`) | `P[i].GravAccel[]` |
