#!/bin/bash
#SBATCH --job-name=AB_BUILD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0:30:00
#SBATCH --mem=8g
# Build batched + per-cell executables. Submitted by submit_ab.sh with --export=ALL,REPO=...,TAG=...,EXTRA="..."
set -e
source "$REPO/examples/agora_disc_star_formation_3d/ab_tests/slurm/ozstar_env.sh"
cd "$REPO"
# EXTRA is a space-separated list of CONFIG_LINEs for build_ab.sh (may be empty)
examples/agora_disc_star_formation_3d/ab_tests/build_ab.sh "$TAG" $EXTRA
