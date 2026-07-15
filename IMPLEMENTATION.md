# SIDM Implementation Reference

Companion to `SIDM.md`. That file explains architecture decisions and *why*; this file is a precise map of *what exists and where* — every struct, function, and call site, so anyone can navigate the codebase without re-reading narrative. Not a substitute for reading the actual header comments in the source, which go into more depth on the trickier pieces.

---

## `src/sidm/sidm.h`

Module header. Included by every SIDM-touching file, plus several files outside `src/sidm/` (see cross-references below).

**`struct DM_Particle_Data`** — the per-DM-particle side-array struct:
```c
typedef struct DM_Particle_Data
{
  MyIDType    PIndex;              // back-reference into P[]
  MyFloat     SidmDensity;
  MyFloat     SidmHsml;
  MyFloat     SidmVelDisp;
  int         SidmNumNgb;
  integertime SidmLastScatterTime;
  int         SidmScatterFlag;
} DM_Particle_Data;
```
`PIndex` mirrors `Star_Particle_Data.PID` mechanically but is named for clarity (it's an array index, not a persistent ID, despite the star-module convention calling the equivalent field `PID`).

**Globals**: `extern DM_Particle_Data *DMSP;`, `extern int NumDM;` — actual (non-`extern`) definitions live in `sidm_density.c`.

**Macros**:
- `DMPS(i)` = `DMSP[P[i].SIDMID]` — given a `P[]` index, get the particle's `DM_Particle_Data` (forward). Mirrors `SPP(i)`.
- `PDMS(i)` = `P[DMSP[i].PIndex]` — given a `DMSP[]` index, get the particle's `particle_data` (backward). Mirrors `PPS(i)`.

**`SIDM_TIMESTEP_SAFETY_FACTOR`** — `#define`d constant (`0.1`), bounds per-step scattering probability in the timestep criterion. Not a `param.txt` entry (see TODO in SIDM).

**Prototypes**: `sidm_density(void)`, `sidm_scatter(int timebin)`. Deliberately not duplicated in `proto.h` — matches the codebase's own convention that module-owned prototypes live in the module's header, not the shared one.

---

## `src/sidm/sidm_tree.h` / `sidm_tree.c`

Dedicated DM-only (`Type==1`) neighbour-search tree, independent of the gravity tree and of `ngbtree` (gas-only). See SIDM §3.1 for why this exists.

**`struct SidmNODE`** (`sidm_tree.h`) — tree node, stripped of mesh-vertex/hydro/RT fields that don't apply to static DM particles. Carries `range_min[3]`/`range_max[3]` (AABB bounds, not gravity's center+len cube), `sibling`/`nextnode`, `father`, `Ti_Current`.

**Globals**: `SidmTree_Nodes`, `SidmTree_MaxPart`, `SidmTree_MaxNodes`, `SidmTree_NumNodes`, `SidmTree_NextFreeNode`, `SidmTree_FirstNonTopLevelNode`, `SidmTree_Nextnode`, `SidmTree_Father`, `SidmTree_Marker`, `SidmTree_DomainNodeIndex`.

**Functions** (`sidm_tree.c`):
- `sidm_treeallocate(void)` — initial allocation, sized off `All.MaxPart` (no `MaxPartDM` exists).
- `sidm_treefree(void)` — silently no-ops if not currently allocated (matches `ngb_treefree()`'s real behaviour; an earlier version incorrectly `terminate()`d here, crashing on `run()`'s first pass — fixed).
- `sidm_treebuild(void)` — driver; construction, mirrors `ngb_treebuild()`.
- `sidm_treebuild_construct(void)` (static) — Peano/Morton-key particle insertion, filtered to `Type==1`, iterating the *full* `NumPart` range (not a contiguous `0..npart` range like gas, since DM isn't contiguous in `P[]`) and skipping non-DM inline.
- `sidm_create_empty_nodes(...)` (static) — builds the empty top-level domain-mirroring skeleton; fully generic, no DM-specific logic, ported near-verbatim from `ngb_create_empty_nodes`.
- `sidm_record_topnode_siblings(...)` (static) — generic, ported near-verbatim.
- `sidm_update_node_recursive(...)` (static) — recursive `range_min`/`range_max` bounding-box propagation up the tree.
- `sidm_exchange_topleafdata(void)` (static) — `MPI_Allgatherv` of each task's local top-leaf ranges to all tasks.
- `sidm_treefind_export_node_threads(int no, int target, int thread_id)` — pseudo-particle export, mirrors `ngb_treefind_export_node_threads`. Reused directly by both `sidm_density.c` and `sidm_scatter.c`.

**Call sites for allocate/free/build**: `run.c`, tied to domain decomposition (SIDM §3.2).

---

## `src/sidm/sidm_density.c`

Hsml-iteration density and velocity-dispersion estimator.

**Globals defined here** (not just declared): `int NumDM;`, `DM_Particle_Data *DMSP;` — the actual definitions for the externs declared in `sidm.h`.

**`typedef struct { MyDouble Pos[3]; MyFloat Hsml; int Firstnode; } data_in;`** — per-particle query, standard `generic_comm_pattern` shape.

**`typedef struct { MyFloat Density; MyFloat VxSum, VySum, VzSum; MyFloat V2Sum; MyFloat Ngb; } data_out;`** — carries *raw* velocity moments, not a pre-computed `VelDisp`. This is the mean-subtraction fix (SIDM §3.8): raw sums are additive across local+imported contributions; a `sqrt()`'d value is not.

**Persistent scratch arrays** (allocated/freed each `sidm_density()` call, in LIFO order): `SidmNgbs`, `Left`, `Right` (Hsml bisection), `SidmVxSum`, `SidmVySum`, `SidmVzSum`, `SidmV2Sum` (raw moment accumulation).

**`particle2in(data_in *in, int i, int firstnode)`** — packs a query.

**`out2particle(data_out *out, int i, int mode)`** — `MODE_LOCAL_PARTICLES` sets, `MODE_IMPORTED_PARTICLES` accumulates (`+=`) onto the scratch arrays and `P[i].SidmDensity`/`SidmNumNgb` directly. Does **not** write `SidmVelDisp` — that's computed once, later, after the bisection loop converges.

**`sidm_density_evaluate(int target, int mode, int threadid)`** (static) — the tree walk. AABB overlap test (not gravity's MAC), cubic-spline kernel (`KERNEL_COEFF_1/2/5`, matching `subfind_density.c`'s convention), accumulates `density_sum`/`vx_sum`/`vy_sum`/`vz_sum`/`v2_sum`/`numngb`. Follows pseudo-particle branches via `sidm_treefind_export_node_threads` (from `sidm_tree.c`).

**`kernel_local(void)` / `kernel_imported(void)`** (static) — standard `generic_comm_pattern` kernels.

**`sidm_density(void)`** — the driver:
1. Left/Right bisection loop (`generic_comm_pattern` each iteration) until `SidmNumNgb` converges to `All.SidmDesNumNgb ± All.SidmDesNumNgbDev`. (`sidm_treebuild()` is **not** called from here — the caller, `run.c`, calls it immediately before, tied to domain decomposition; see SIDM §3.2.)
2. **After** the loop converges: a dedicated pass computing `P[i].SidmVelDisp = sqrt(max(0, V2Sum/Density - (VxSum²+VySum²+VzSum²)/Density²) / 3)` — the actual mean-subtracted dispersion, using the *final* converged pass's accumulated moments.
3. MPI-reduced diagnostic print (min/max/mean of `Density`/`Hsml`/`NumNgb`) — still the primary way to inspect these values until `A_DMSP` I/O was wired (now redundant with snapshot output, but left in place).

**Call site**: `run.c`, immediately after `sidm_treebuild()`, both tied to domain decomposition.

---

## `src/sidm/sidm_scatter.c`

Monte Carlo pairwise elastic scattering/kick, local and cross-task partners.

**`enum { MAX_CAND_TOTAL = 128, MAX_REMOTE_RESPONSE = 16 };`** — fixed caps (SIDM TODO: not dynamically sized).

**`struct sidm_candidate`** — unified local/remote candidate: `is_remote`, `index` (local `P[]` index, or remote `P[]` index on `remote_task` if `is_remote`), `remote_task`, `r`, `p_ij`, `vel[3]` (only populated/used if remote).

**`SidmCandBuf`** (flattened `[NumPart][MAX_CAND_TOTAL]`) / **`SidmCandCount`** — per-particle candidate list and count, built during the comm-pattern walk, consumed by the accept/reject pass afterward.

**`sidm_dt_code_units(int timebin)`** (static) — physical Δt for the given timebin, via `get_time_difference_in_Gyr()` (SIDM §3.7 — **not** a raw `Timebase_interval` multiply, that would be `Δ(ln a)` in a comoving run, not a duration).

**`sidm_random_unit_vector(double e[3])`** (static) — isotropic kick direction.

**`sidm_check_conservation(...)`** (static) — shared momentum/energy diagnostic, called by both kick paths below. Prints a warning if relative error exceeds `1e-10`; tracks running max for the summary print.

**`sidm_apply_kick_local(int i, int j)`** (static) — both particles local. Direct port of the pre-cross-task-support version: CM-frame isotropic kick, Eq. 20.

**`struct sidm_kick_poke`** — `{ int remote_index; MyFloat new_vel[3]; integertime scatter_time; }` — the one-way "set this velocity" instruction for a remote kick.

**`sidm_apply_kick_cross_task(int i, const sidm_candidate *cand, sidm_kick_poke **poke_buf, int *poke_count, int poke_capacity)`** (static) — computes both new velocities using `cand->vel` (already fetched during gathering), applies the local half immediately, queues the remote half into `poke_buf[cand->remote_task]`. `terminate()`s if `remote_task` is out of range (defensive — should be unreachable if self-reporting worked) or if the per-task poke buffer overflows its provisioned capacity.

**`data_in`** (candidate-gathering query) — `Pos[3]`, `Hsml`, `SelfID` (needed for the remote side's lower-ID-initiates filter), `Firstnode`.

**`data_out`** (candidate-gathering response) — `ncand`, `remote_task` (**self-reported** `ThisTask` — the comm framework doesn't otherwise tell the originating side which task a response came from), then per-candidate `id[]`, `vel[][3]`, `r[]`, `remote_index[]`, up to `MAX_REMOTE_RESPONSE`.

**`particle2in`** — packs a query, mirrors `sidm_density.c`.

**`out2particle(data_out *out, int i, int mode)`** — `MODE_LOCAL_PARTICLES` is a no-op (local candidates are appended directly into `SidmCandBuf` *during* the walk, not here). `MODE_IMPORTED_PARTICLES` computes `p_ij` for each returned candidate (using `i`'s own `Mass`/`Vel`, available on the originating task) and appends into `SidmCandBuf[i]`.

**`sidm_scatter_walk(int target, int mode, int threadid)`** (static) — the shared walk, used both for a particle's own local search (`MODE_LOCAL_PARTICLES`, appending directly into `SidmCandBuf`) and for serving a remote query (`MODE_IMPORTED_PARTICLES`, packing into the `data_out` response, self-reporting `remote_task = ThisTask`). Unlike `sidm_density_evaluate`, this one **genuinely follows pseudo-particle branches** via `sidm_treefind_export_node_threads` rather than skipping them — this is the actual cross-task-support change.

**`kernel_local(void)` / `kernel_imported(void)`** (static) — standard shape; `kernel_local` additionally filters to `P[i].TimeBinGrav == sidm_scatter_current_timebin` (see below) and resets `SidmCandCount[i] = 0` before walking.

**`sidm_scatter(int timebin)`** — the driver:
1. Computes `dt_i`, stores it in file-static `sidm_scatter_current_dt_i` (read by `out2particle`/the walk).
2. Allocates `SidmCandBuf`/`SidmCandCount`, runs `generic_comm_pattern(kernel_local, kernel_imported)` — candidate gathering, local + remote.
3. Allocates per-destination-task pending-poke buffers (`poke_buf[t]`, fixed capacity `TimeBinsGravity.NActiveParticles` per task).
4. Iterates active particles at exactly `timebin` (matching `P[i].TimeBinGrav == timebin` — **not** the cumulative hierarchical active list gravity's own kick uses, which would re-roll a particle's scatter dice once per coarser timebin pass): accept/reject + partner selection over the combined candidate list; dispatches to `sidm_apply_kick_local` or `sidm_apply_kick_cross_task`.
5. **Kick delivery**: a plain per-task `MPI_Sendrecv` loop (count exchange on `TAG_N`, payload on `TAG_DMDATA`) delivering queued pokes; the receiving side applies `P[ri].Vel = new_vel`, `DMPS(ri).SidmLastScatterTime`/`SidmScatterFlag`.
6. Frees everything in LIFO order; prints the MPI-reduced `SIDM_SCATTER:` diagnostic line, now including a separate `cross_task=N` count.

**Call site**: `do_gravity_hydro.c`, immediately after gravity's own per-timebin kick loop, inside the `HIERARCHICAL_GRAVITY` branch only (see SIDM TODO for the non-hierarchical gap).

---

## Files touched outside `src/sidm/`

### `src/main/allvars.h`
- `particle_data.SIDMID` (`MyIDType`) — forward reference into `DMSP[]`.
- `global_data_all_processes` (`All.`): `SidmDesNumNgb`, `SidmDesNumNgbDev`, `SidmCrossSection`.
- `enum e_typelist`: `DM_ONLY = 2`.
- `enum arrays`: `A_DMSP`.
- `enum iofields`: `IO_SIDM_DENSITY`, `IO_SIDM_HSML`, `IO_SIDM_NUMNGB`, `IO_SIDM_VELDISP`.

### `src/main/run.c`
- Includes `sidm.h`.
- `calculate_non_standard_physics_with_valid_gravity_tree(void)` — back to an empty stub; no longer the SIDM hook (see SIDM §3.2).
- Two call sites (both domain-decomposition-gated blocks): `sidm_treeallocate()`/`sidm_treefree()` alongside `ngb_treeallocate()`/`ngb_treefree()`; `sidm_treebuild()` + `sidm_density()` immediately after `ngb_treebuild(NumGas)`.

### `src/io/io.c`
- `init_field()` — `A_DMSP` offset computation (`(size_t)pointer_to_field - (size_t)DMSP`).
- `fill_write_buffer()` — full `A_DMSP` write-path block, mirroring the `A_S` block exactly (uses `PDMS(pindex).Type` instead of `PPS(pindex).Type`).

### `src/io/io_fields.c`
- Four `init_field`/`init_units`/`init_snapshot_type` triples for `IO_SIDM_DENSITY`/`HSML`/`NUMNGB`/`VELDISP`, now pointing at `&DMSP[0].Sidm...` with `A_DMSP` (previously pointed at `&P[0].Sidm...` with `A_P`, from before the `DMSP[]` migration).

### `src/io/read_ic.c`
- `empty_read_buffer()` — `A_DMSP` case in the array-source switch (`array_pos = DMSP + n`).
- Two per-type particle-counting loops (there are two nearly-identical ones in this file) — `NumDM += n_for_this_task` when `type==1`.
- `NumDM = 0` reset alongside `NumStars`/`NumBhs`.
- The actual index-assignment loop (mirrors `STARS`/`BLACKHOLES`): `P[id].SIDMID = jd; DMSP[jd].PIndex = id;` for every `Type==1` particle, in ID order. **This was the fix for the bus-error crash** — nothing established these indices for a freshly-read IC before this was added.

### `src/io/restart.c`
- `in(&NumDM, modus)` + `byten(&DMSP[0], NumDM*sizeof(DM_Particle_Data), modus)`, positioned between the `STAR_FEEDBACK_ACTIVE` block and `BLACKHOLES`, mirroring `SP[]`'s own checkpoint block exactly.

### `src/io/parameters.c`
- `param.txt` registration for `SidmDesNumNgb`, `SidmDesNumNgbDev`, `SidmCrossSection` (all `REAL` type), mirroring the `StarDesNgb`/`StarDesDev` pattern.

### `src/utils/allocate.c`
- `DMSP` initial allocation (`mymalloc_movable`, sized `All.MaxPart`), zero-init (`memset`), and reallocation inside `reallocate_memory_maxpart()` (lockstep with `P[]`'s own growth — no separate `reallocate_memory_maxpartdm()`-style function exists or is needed, since `DMSP[]` doesn't have its own independently-changing size parameter the way `SP[]`/`BhP[]` do).

### `src/utils/tags.h`
- `TAG_DMDATA` (value `17`) — used consistently on both sides of every SIDM `Sendrecv`/`Isend`/`Irecv` call (deliberately avoiding the mismatched-tag pattern found, but not fixed, in the pre-existing `STARS` code — see SIDM §6 TODO).

### `src/domain/domain.h`, `domain_vars.c`, `domain_counttogo.c`
- `toGoDM`/`toGetDM` — declared (`domain.h`), defined and allocated/freed in correct LIFO order (`domain_vars.c`), counted per-particle and exchanged via `MPI_Alltoall` (`domain_counttogo.c`) — full mirror of `toGoStars`/`toGetStars`.

### `src/domain/domain_exchange.c`
The largest single mirror of the `STARS` pattern — roughly 20 distinct touch points in one function:
- Variable declarations (`count_dm`, `offset_dm`, `count_recv_dm`, `offset_recv_dm`, `dmBuf`).
- Allocation of the above, in the same relative order as `STARS`'s own (for correct LIFO free later).
- Initial offset computation, count totals.
- `dmBuf` allocation (`count_togo_dm * sizeof(DM_Particle_Data)`).
- Per-particle export-packing branch (`else if(P[n].Type == 1)`), between the `STARS` and `BLACKHOLES` branches.
- Local swap-removal: the new `Type==1` branch itself, **plus** added `DMPS`/`PIndex` cross-checks to the *existing* `GAS`/`STARS`/`BLACKHOLES` branches (they didn't check for a shuffled-in DM particle before this).
- The `count_totget` index-shift block (`DMSP[i].PIndex += count_totget`) — accounts for gas particles being inserted at the front of `P[]`.
- `count_recv_dm`/`offset_recv_dm` computation.
- The actual exchange: `Sendrecv` (using one consistent tag on both sides — **not** replicating the pre-existing `TAG_BHDATA`/`TAG_STARDATA` mismatch found in the `STARS` code's own equivalent call, see SIDM §6), `Irecv`, `Isend`, and the `myMPI_Alltoallv` alternate-path variant.
- Import-side index assignment (`P[i].SIDMID = j; DMSP[j].PIndex = i;`) for newly-arrived DM particles.
- All frees, in correct LIFO order (`dmBuf` between `bhBuf` and `sBuf`; the four `count`/`offset` arrays between `BLACKHOLES`'s and `STARS`'s own).

### `src/domain/peano.c`
- DM back-reference fixup (`DMPS(dest).PIndex = dest`) inside the Peano-key cycle-sort, alongside the existing `STARS`/`BLACKHOLES` checks. Runs on essentially every domain decomposition — a real, frequently-hit gap before this was added.

### `src/domain/domain_rearrange.c`
- DM back-reference fixup in the gas-cell-elimination branch (when the last particle shuffles into a freed slot).

### `src/mesh/refinement.c`
- `move_collisionless_particle()` — DM back-reference fixup, alongside `STARS`/`BLACKHOLES`. This function's own doc comment confirms it moves *any* collisionless particle when gas refinement shifts the gas/non-gas boundary — includes DM by definition, and this can happen more often than domain decomposition itself in a real AGORA run.

### `src/time_integration/timestep.c`
- `get_timestep_gravity(int p)` — SIDM criterion added right after the `EXTERNALGRAVITY` block: `dt_sidm = SIDM_TIMESTEP_SAFETY_FACTOR / (rho_phys · σ/m · v_rel)`, `rho_phys` via `All.cf_a3inv`, `v_rel` via `All.cf_atime` (matching this file's own existing precedent in `get_timestep_hydro`, not independently derived). Skipped if `DMPS(p).SidmDensity <= 0` (density never computed for this particle yet).

### `src/time_integration/do_gravity_hydro.c`
- `sidm_scatter(timebin)` call, right after the per-timebin `kick_particle` loop, inside the `HIERARCHICAL_GRAVITY` branch. The non-hierarchical branch is explicitly *not* wired (see SIDM TODO) — flagged with a comment rather than silently left blank.

### `Config.sh` / `Makefile`
- `SIDM` flag; `sidm/sidm_density.o`, `sidm/sidm_tree.o`, `sidm/sidm_scatter.o` build rules; `sidm/sidm.h`, `sidm/sidm_tree.h` in `INCL`.

---

## Quick lookup: "where is X actually computed?"

| Quantity | Computed in | Persisted in |
|---|---|---|
| `SidmDensity` | `sidm_density_evaluate()` (kernel-weighted sum) | `DMPS(i).SidmDensity` |
| `SidmHsml` | `sidm_density()`'s bisection loop | `DMPS(i).SidmHsml` |
| `SidmVelDisp` | `sidm_density()`, **after** the bisection loop converges (mean-subtracted) | `DMPS(i).SidmVelDisp` |
| `SidmNumNgb` | `sidm_density_evaluate()` | `DMPS(i).SidmNumNgb` |
| Scattering probability `P_ij` | `sidm_scatter_walk()` (local) / `out2particle()` (remote) | transient, in `SidmCandBuf` only |
| `SidmLastScatterTime`/`SidmScatterFlag` | `sidm_apply_kick_local()` / `sidm_apply_kick_cross_task()` / kick-delivery receive loop | `DMPS(i).Sidm*` |
| `σ/m` in code units | *you*, from `param.txt` `SidmCrossSection` | `All.SidmCrossSection` — see SIDM §3.4 for the conversion |
