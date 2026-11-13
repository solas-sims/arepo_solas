#!/bin/bash
# Arepo-Solas MPI run on Setonix (Pawsey), CPU "work" partition.
# DRAFT, not yet run: check ACCOUNT, the modules and the library paths before use.
#
# Usage (from the run directory, which holds param.txt, ICs/ and output_list.txt; logs/ is created):
#   sbatch --account=<project> /path/to/submit_setonix.sbatch
#   sbatch --account=<project> --export=ALL,EXEC=$PWD/Arepo-CosmoBox,NTASKS=64 submit_setonix.sbatch
#
# Restarts: the script starts from the ICs if output/restartfiles/restart.0 does not exist, otherwise it
# restarts from the restart files (RestartFlag 1). Arepo stops cleanly and writes restart files once it has used
# 85% of TimeLimitCPU (param.txt), so TimeLimitCPU must stay below ~1.17 x the wall time requested here.
# To continue a long run, chain jobs:  sbatch --dependency=afterany:<jobid> submit_setonix.sbatch
#
# Memory: Setonix work nodes have 128 cores and roughly 230 GB usable, about 1.8 GB per core. param.txt has
# MaxMemSize = 2500 (MB per MPI task), so do not run 128 tasks per node with that value: use 64 tasks per node
# (the default below) or lower MaxMemSize.

#SBATCH --job-name=arepo_solas
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=logs/arepo_%j.out
#SBATCH --error=logs/arepo_%j.err
# The account is not set here; pass --account=<project> on the sbatch line, or uncomment and edit:
##SBATCH --account=<project>

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

# ---- settings (override with --export=ALL,NAME=value) ------------------------------------------------------
EXEC="${EXEC:-$PWD/Arepo}"            # path to the compiled executable
PARAM="${PARAM:-param.txt}"
OUTDIR="${OUTDIR:-./output}"          # must match OutputDir in the parameter file
SNAPNUM="${SNAPNUM:-}"                # set to restart from a snapshot (RestartFlag 2) instead of restart files

# ---- environment: use the same modules as at build time ----------------------------------------------------
# The batch job inherits the submitting shell's environment (PrgEnv-gnu, cray-mpich, xpmem as at build time);
# the libraries are loaded explicitly so the job does not depend on what was loaded when you submitted.
set +u   # the module system reads unset variables
module load gsl/2.8 fftw/3.3.10 hdf5/1.14.5-parallel-api-v112
set -u
# and keep the build-time library paths:
export LD_LIBRARY_PATH="$HOME/software/grackle_solas/lib:$HOME/software/gmp-6.3.0/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1

: "${PAWSEY_HDF5_HOME:?PAWSEY_HDF5_HOME is not set: load the same hdf5/gsl/fftw modules as for the build}"
[ -x "$EXEC" ] || { echo "executable not found or not executable: $EXEC" >&2; exit 1; }
[ -f "$PARAM" ] || { echo "parameter file not found: $PARAM" >&2; exit 1; }

# ---- start or restart --------------------------------------------------------------------------------------
if [ -n "$SNAPNUM" ]; then
  FLAG="2 $SNAPNUM"
elif [ -e "$OUTDIR/restartfiles/restart.0" ]; then
  FLAG="1"
else
  FLAG=""
fi

mkdir -p "$OUTDIR"
echo "Arepo run: ${SLURM_NTASKS} MPI tasks on ${SLURM_JOB_NUM_NODES} node(s), restart flag: '${FLAG:-0 (fresh start)}'"
echo "exec=$EXEC param=$PARAM host=$(hostname) started $(date)"
module list 2>&1 || true

# word splitting of $FLAG is intended (it may be "2 <snapshot>")
STATUS=0
srun -n "$SLURM_NTASKS" "$EXEC" "$PARAM" $FLAG || STATUS=$?

echo "Arepo exited with status ${STATUS} at $(date)"
exit $STATUS
