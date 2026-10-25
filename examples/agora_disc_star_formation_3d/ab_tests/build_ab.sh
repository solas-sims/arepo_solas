#!/bin/bash
# Build the batched and per-cell (GRACKLE_PER_CELL) executables from the same source, for A/B tests.
#
# usage (from the repository root):
#   examples/agora_disc_star_formation_3d/ab_tests/build_ab.sh <tag> [CONFIG_LINE ...]
#
# Extra CONFIG_LINEs are appended to the example's Config.sh (a line with the same KEY replaces the existing one), e.g.
#   build_ab.sh chem1 GRACKLE_CHEMISTRY=1
#   build_ab.sh limiter COOLING_TIMESTEP_LIMITER
#   build_ab.sh jeans JEANS_SF JEANS_MASS_BASED   # also remove AGORA_SF: KEY= (empty value) deletes a line, e.g. AGORA_SF=
# Output: <example>/ab/<tag>/Arepo_batched and Arepo_percell
set -e
TAG=$1
shift
EX=examples/agora_disc_star_formation_3d
OUT=$EX/ab/$TAG
mkdir -p "$OUT"

for MODE in batched percell; do
  CFG=$OUT/Config_$MODE.sh
  cp $EX/Config.sh "$CFG"
  for L in "$@"; do
    KEY=${L%%=*}
    VAL=${L#*=}
    sed -i.bak -E "/^${KEY}([[:space:]=]|$)/d" "$CFG" && rm -f "$CFG.bak"
    [ "$L" = "$KEY=" ] || echo "$L" >> "$CFG"
  done
  [ $MODE = percell ] && echo "GRACKLE_PER_CELL" >> "$CFG"
  make -j8 CONFIG="$CFG" BUILD_DIR="$OUT/build_$MODE" EXEC="$OUT/Arepo_$MODE"
done
