#!/bin/bash
#SBATCH --job-name=AB_RUN
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=8g
# Run one executable. --ntasks is set at submission (sbatch --ntasks=N). Needs --export=ALL,REPO=...,TAG=...,MODE=...,TMAX=...[,TSNAP=...]
set -e
source "$REPO/examples/agora_disc_star_formation_3d/ab_tests/slurm/ozstar_env.sh"
cd "$REPO/examples/agora_disc_star_formation_3d"
LAUNCHER="srun -n" ab_tests/run_ab.sh "$TAG" "$MODE" "$SLURM_NTASKS" "$TMAX" $TSNAP
