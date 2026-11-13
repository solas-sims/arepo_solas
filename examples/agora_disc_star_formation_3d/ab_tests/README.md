# A/B tests: batched vs per-cell Grackle

Checks that the array-based Grackle path (`GrackleCoolCells()` / `GrackleCoolingTimes()`) gives the same answer as the old
one-cell-at-a-time path, and measures the speed-up. Both executables are built from the same source; `GRACKLE_PER_CELL`
switches `cool_cells()` (and the eEOS cooling-time pass) back to the per-cell calls.

Run everything from the example directory, after `./get_ics.sh` and `./get_tables.sh`. `build_ab.sh` is run from the repository root.

```bash
# 1. build (repository root)
examples/agora_disc_star_formation_3d/ab_tests/build_ab.sh base

# 2. run (example directory); LAUNCHER defaults to "mpirun -np", use e.g. LAUNCHER="srun -n" on Setonix/OzSTAR
cd examples/agora_disc_star_formation_3d
ab_tests/run_ab.sh base batched 1 0.01
ab_tests/run_ab.sh base percell 1 0.01

# 3. compare (needs python with h5py + numpy)
ab_tests/compare_snapshots.py ab/base/run_batched_n1/output ab/base/run_percell_n1/output
```

`TimeMax` is in code time (1 unit = 0.978 Gyr), so 0.01 is about 10 Myr; the full run is 0.5.

## Test matrix

| # | What | How | Expect |
|---|---|---|---|
| 1 | First step (`dt = 0`) | `run_ab.sh base <mode> 1 0.0005 0.0001` | identical |
| 2 | Short run, 1 and 4 tasks | `run_ab.sh base <mode> {1,4} 0.01` | identical (`solve_chemistry` is per-cell independent; only OpenMP in Grackle could introduce round-off, then use `--rtol 1e-10`) |
| 3 | Variants | `build_ab.sh chem1 GRACKLE_CHEMISTRY=1`; `build_ab.sh limiter COOLING_TIMESTEP_LIMITER`; `build_ab.sh jeans AGORA_SF= JEANS_SF JEANS_MASS_BASED`; `build_ab.sh eeos AGORA_SF= EEOS_SF` | identical |
| 4 | Full run to TimeMax 0.5 | `run_ab.sh base <mode> N 0.5 0.02`, then `agora_disc_diagnostics.py` and the SFR history (`output/sfr.txt`) on both | agree statistically; late-time divergence is expected if test 2 is not bit-exact |
| 5 | Performance | compare the `sfrcool` row of the final cumulative table in `output/cpu.txt` for both modes at 1, 4, 8+ tasks | batched faster |
| 6 | Robustness | restart from a mid-run restart file; build with `-fsanitize=address` | no crashes, no leaks |

Notes:

- Test 3 `eeos` is the case least exercised: it is eEOS + Grackle, which also uses the cooling-time sign convention
  (Grackle returns negative cooling times). Check `Sfr` and `InternalEnergy` carefully there.
- Particles are matched by `ParticleIDs`. If refinement/derefinement or star formation diverges between runs the particle
  sets differ and the comparison reports it rather than comparing mismatched cells; that is itself the signal of a difference.
- `ab/` is git-ignored.

## SLURM (OzSTAR/NT)

`slurm/` holds batch scripts built from the `AGORA_LOWRES_TEST` job script. Edit `slurm/ozstar_env.sh` first: it loads
`gcc/11.2.0 openmpi/4.1.1 hdf5/1.12.1 fftw/3.3.10 gsl/2.7`, and you must add a python with numpy + h5py (for the comparison) and make
sure the Makefile can find Grackle (`GRACKLE_DIR` in `Makefile.systype`).

```bash
cd examples/agora_disc_star_formation_3d && ./get_ics.sh && ./get_tables.sh   # on a node with internet access
ab_tests/slurm/submit_ab.sh base 4 0.01 0.0025                     # build -> batched + per-cell runs -> compare
ab_tests/slurm/submit_ab.sh chem1 4 0.01 0.0025 -- GRACKLE_CHEMISTRY=1
SKIP_BUILD=1 ab_tests/slurm/submit_ab.sh base 8 0.01 0.0025        # reuse executables, new task count
RTOL=1e-10 ab_tests/slurm/submit_ab.sh base 4 0.01                 # allow round-off
```

`submit_ab.sh` chains the jobs with `--dependency=afterok`. Logs are in `ab/<tag>/slurm_logs/`, the comparison result in
`ab/<tag>/compare_n<ntasks>.log`. Config overrides after `--` must not contain commas. Run times: `run_ab.sbatch` asks for 4 h
and 8 GB per task; edit the `#SBATCH` lines for the full 0.5 run.
