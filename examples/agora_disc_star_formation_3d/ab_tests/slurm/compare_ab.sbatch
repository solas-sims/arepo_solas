#!/bin/bash
#SBATCH --job-name=AB_COMPARE
#SBATCH --ntasks=1
#SBATCH --time=0:30:00
#SBATCH --mem=16g
# Compare batched vs per-cell output. Needs --export=ALL,REPO=...,TAG=...,NTASKS=...[,RTOL=...]
source "$REPO/examples/agora_disc_star_formation_3d/ab_tests/slurm/ozstar_env.sh"
cd "$REPO/examples/agora_disc_star_formation_3d/ab"/"$TAG"
"$REPO/examples/agora_disc_star_formation_3d/ab_tests/compare_snapshots.py" \
  run_batched_n"$NTASKS"/output run_percell_n"$NTASKS"/output --rtol "${RTOL:-0}" | tee compare_n"$NTASKS".log
exit "${PIPESTATUS[0]}"
