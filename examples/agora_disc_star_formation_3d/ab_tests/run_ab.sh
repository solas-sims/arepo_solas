#!/bin/bash
# Run one A/B executable on the AGORA disc for a short time.
#
# usage (from the example directory, after get_ics.sh and get_tables.sh):
#   ab_tests/run_ab.sh <tag> <batched|percell> <ntasks> <TimeMax> [TimeBetSnapshot]
# PARAM_OVERRIDES: space-separated Key=value pairs replacing/adding entries in param.txt, e.g. PARAM_OVERRIDES="MinVolume=5e-4".
# Launcher: set LAUNCHER, e.g. LAUNCHER="srun -n" (Setonix/OzSTAR) -- default "mpirun -np".
# Output: ab/<tag>/run_<mode>_n<ntasks>/output
# TimeMax is in code time (UnitTime = 0.978 Gyr, so 0.01 ~ 9.8 Myr; the full run is 0.5).
set -e
TAG=$1; MODE=$2; N=$3; TMAX=$4; TSNAP=${5:-$(python3 -c "print($TMAX/4)")}
HERE=$(pwd)
RUN=$HERE/ab/$TAG/run_${MODE}_n${N}
mkdir -p "$RUN/output"
cd "$RUN"
for D in ICs grackle_data star_tables; do ln -sfn "$HERE/$D" "$D"; done
sed -E "s|^TimeMax .*|TimeMax $TMAX|; s|^TimeBetSnapshot .*|TimeBetSnapshot $TSNAP|; s|^TreecoolFile .*|TreecoolFile $HERE/../../TREECOOL_ep|" \
  "$HERE/param.txt" > param.txt
for KV in $PARAM_OVERRIDES; do
  K=${KV%%=*}; V=${KV#*=}
  if grep -qE "^$K[[:space:]]" param.txt; then
    sed -i.bak -E "s|^$K[[:space:]].*|$K $V|" param.txt && rm -f param.txt.bak
  else
    echo "$K $V" >> param.txt
  fi
done
${LAUNCHER:-mpirun -np} "$N" "$HERE/ab/$TAG/Arepo_$MODE" param.txt 2>&1 | tee run.log
