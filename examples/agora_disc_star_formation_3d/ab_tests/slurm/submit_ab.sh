#!/bin/bash
# Submit build -> (batched, per-cell) runs -> comparison as a dependency chain on OzSTAR/NT.
#
# usage: submit_ab.sh <tag> <ntasks> <TimeMax> [TimeBetSnapshot] [-- CONFIG_LINE ...]
#   submit_ab.sh base 1 0.01
#   submit_ab.sh base 4 0.01 0.0025
#   submit_ab.sh chem1 4 0.01 0.0025 -- GRACKLE_CHEMISTRY=1
# Set SKIP_BUILD=1 to reuse executables from an earlier build with the same tag; RTOL=1e-10 to allow round-off.
# The example needs ICs/ and tables (get_ics.sh, get_tables.sh) in examples/agora_disc_star_formation_3d first.
set -e
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
TAG=$1; NTASKS=$2; TMAX=$3; shift 3
TSNAP=""
if [ $# -gt 0 ] && [ "$1" != "--" ]; then TSNAP=$1; shift; fi
[ "$1" = "--" ] && shift
EXTRA="$*"

EX=$REPO/examples/agora_disc_star_formation_3d
[ -f "$EX/ICs/agora_lowres_arepo_ic_without_u-with-grid.hdf5" ] || { echo "ICs missing: run get_ics.sh in $EX"; exit 1; }
[ -f "$EX/grackle_data/CloudyData_UVB=FG2011_shielded.h5" ] || { echo "Grackle table missing: run get_tables.sh in $EX"; exit 1; }
mkdir -p "$EX/ab/$TAG/slurm_logs"
LOG="$EX/ab/$TAG/slurm_logs"
EXPORT="ALL,REPO=$REPO,TAG=$TAG,TMAX=$TMAX,TSNAP=$TSNAP,NTASKS=$NTASKS,RTOL=${RTOL:-0}"

DEP=""
if [ -z "$SKIP_BUILD" ]; then
  B=$(sbatch --parsable --export="$EXPORT,EXTRA=$EXTRA" -o "$LOG/build_%j.out" "$HERE/build_ab.sbatch")
  DEP="--dependency=afterok:$B"
  echo "build:   $B"
fi

R1=$(sbatch --parsable $DEP --ntasks="$NTASKS" --export="$EXPORT,MODE=batched" -o "$LOG/run_batched_%j.out" "$HERE/run_ab.sbatch")
R2=$(sbatch --parsable $DEP --ntasks="$NTASKS" --export="$EXPORT,MODE=percell" -o "$LOG/run_percell_%j.out" "$HERE/run_ab.sbatch")
echo "batched: $R1"
echo "percell: $R2"
C=$(sbatch --parsable --dependency=afterok:$R1:$R2 --export="$EXPORT" -o "$LOG/compare_%j.out" "$HERE/compare_ab.sbatch")
echo "compare: $C   (log: $LOG/compare_$C.out, ab/$TAG/compare_n$NTASKS.log)"
