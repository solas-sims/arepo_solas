# `pawsey-pulse-gpu` branch: CI and branch protection

This documents the branch set up on 2026-09-09 as a stable base for the Pawsey PULSE
team's work porting stellar/BH feedback code to GPU, and the two CI workflows added to
gate it.

## Branch

`pawsey-pulse-gpu`, branched from `main` at commit `993afd2` ("Add small cosmological
box example"). Nothing else changed on it besides the CI setup described below — it's
meant as a known-good starting point, not a feature branch of its own.

## `build.yml` — compile matrix

Two jobs:

- `discover-examples` parses `test.sh` for its `TESTS+="..."` lines (skipping commented-out
  ones) and turns the result into a matrix, so the list of examples compiled here can't
  drift out of sync with what `test.sh` itself exercises.
- `compile` builds each example's `Config.sh` with `make` — compile only, no simulation run.
  Runs on `ubuntu-latest` with `libopenmpi-dev`, `libgsl-dev`, `libhdf5-serial-dev`,
  `libfftw3-dev`, `libhwloc-dev` installed via `apt`. Fifteen examples as of this writing.

A third job, `build-summary`, depends on `compile` and reports a single pass/fail result
as the check named **`Build (all examples)`**. This exists because the matrix produces one
check per example (`compile (wave_1d)`, `compile (noh_3d)`, etc.), and that list shifts
every time an example is added or removed — not something you want to have to update in
branch protection settings each time.

Running this caught a real, pre-existing bug: `examples/cosmo_box_star_formation_3d/Config.sh`
was missing a star-formation model flag that `USE_SFR` now requires (one of `EEOS_SF`,
`AGORA_SF`, `JEANS_SF`). This predates that requirement being added and had nothing to do
with GPU work — fixed by adding `EEOS_SF`, matching the identical SF parameters already
used in `galaxy_merger_star_formation_3d`.

Triggers: push and pull request against `pawsey-pulse-gpu` and `main`, plus manual dispatch.

## `agora-smoke.yml` — build-and-run check

Builds and runs `examples/agora_disc_star_formation_3d`, which exercises the
`STARS`/`STAR_PARTICLES`/`WINDS`/`SUPERNOVAE` feedback path the GPU work will target.
Reports as the check named **`agora-smoke`**.

Steps:

1. Install build dependencies (same as above, plus `cmake`).
2. Build Grackle from `solas-sims/grackle`, pinned to commit `7fbc9bb` (that fork has no
   tags yet). Uses CMake, not the classic `./configure` build — the latter is deprecated
   upstream as of Grackle 3.4 and scheduled for removal in 3.5. Cached across runs.
3. Download the AGORA disc ICs and star-feedback tables via the example's own
   `get_ics.sh`/`get_tables.sh` (Dropbox links). Cached across runs.
4. Build Arepo against `examples/agora_disc_star_formation_3d/Config.sh`.
5. Run the example for a *target number of sync-points* rather than a fixed run time.

### Why sync-points, not a timer

Two earlier approaches to bounding this run didn't work, and both are still explained
inline in the workflow file for anyone tempted to redo them:

- **A small `TimeMax`** derived from `target_steps * MaxSizeTimestep`. Doesn't work because
  hierarchical/individual timestepping means the number of sync-points per unit of
  simulated time isn't predictable up front. One run raced through 700+ sync-points in 43
  minutes without getting anywhere near the reduced `TimeMax`, and had to be cancelled.
- **A fixed `TimeLimitCPU`** (tried 60s). Doesn't work because startup — reading the ICs,
  domain decomposition, building the mesh and gravity tree — eats a large and apparently
  variable share of any fixed wall-clock budget. One run only reached 2 sync-points in 60s.

What actually works: launch Arepo in the background, poll its log for the real sync-point
count, and once the target is hit, drop an `output/stop` file. Arepo checks for that file
at every sync-point (`check_for_interruption_of_run()` in `src/main/run.c`) and exits
cleanly (exit code 0), writing a restart file on the way out. This was checked against a
stand-in shell script that mimics the log format and stop-file check before trusting it
against the real binary.

Defaults: 15 sync-points on push/PR, 100 sync-points on the weekly schedule run.
Configurable via `workflow_dispatch` (`target_sync_points` input). A generous outer wait
(scaled to the target, floored at 600s) guards against the run stalling; the job's overall
90-minute timeout is the final backstop.

### Pass/fail check

After the run: at least 5 sync-points must be logged, and at least one snapshot file must
exist under `output/`. On failure, `output/`, `param_ci.txt`, and `run.log` are uploaded as
a build artifact for inspection.

### Expected noise in the log

Lines like:

```
WARNING: Code warning on task=0, function do_derefinements(), file src/mesh/voronoi/voronoi_derefinement.c, line 594: strange: fac=-4.4e-36
```

are not a bug and not something introduced by this branch. They come from pre-existing
upstream code (commit `60bd71b`, 2022) in `do_derefinements()`: when a cell is de-refined,
its mass/momentum is redistributed to neighbours proportional to a volume fraction; a
degenerate/sliver Voronoi triangle can occasionally produce a tiny negative fraction from
floating-point round-off (the magnitudes seen are noise, on the order of 1e-25 to 1e-36,
not real negative volumes). The code already handles this — it warns, clamps the fraction
to zero, and continues (there's a commented-out `terminate("strange")` in the same block
that was deliberately left disabled). Do not treat this warning as a smoke-test failure
signal.

### Triggers

- Push and pull request against `pawsey-pulse-gpu`.
- Weekly, Monday 03:17 UTC — this job depends on external Dropbox links that aren't under
  our control, so it isn't run on every PR to avoid blocking merges on their availability.
- Manual dispatch, with `target_sync_points` overridable.

## Branch protection on `pawsey-pulse-gpu`

Set via the GitHub API on 2026-09-09:

- Required status checks: `Build (all examples)`, `agora-smoke`.
- `strict: true` — the branch must be up to date with the target branch before these
  checks are considered valid, so a PR can't merge on a stale pass.
- Force-pushes and branch deletion are blocked.
- Not enforced on repo admins (`enforce_admins: false`).

No required PR review count was set — add one later if the team wants it once there's
more than one regular contributor pushing to this branch.
