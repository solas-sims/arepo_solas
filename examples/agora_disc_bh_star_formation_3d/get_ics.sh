#!/bin/bash
# Prepares the AGORA isolated-disc initial conditions *with a central black hole*
# (1e6 Msun PartType5 at the disc centre) in ./ICs/.
#
#   ./get_ics.sh lowres                 downloads the low-res IC and adds the BH
#   ./get_ics.sh medres [PATH_TO_IC]    adds the BH to the medium-res IC you supply
#                                       (default ./ICs/agora_galaxy_med_arepo-with-grid.hdf5)
#
# Needs python3 with h5py and numpy. BH mass: BH_MASS_MSUN (default 1e6).
#
# The medium-res IC is not downloaded: it is the AGORA MED set already passed through ADDBACKGROUNDGRID.
# In both ICs the gas InternalEnergy is set to 0, so Arepo initialises it from InitGasTemp (1e4 K) in param_*.txt.
# The medium-res copy is also changed as follows, to match the low-res IC (see README):
#   - header UnitMass_in_g 1.989e43 -> 1.989e42 (its masses are in 1e9 Msun)
#   - gas InternalEnergy set to 0 (as for low-res)
#   - the BH is placed at BoxSize/2, which is the disc centre in this file

set -euo pipefail

RES="${1:-lowres}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IC_DIR="./ICs"
BH_MASS_MSUN="${BH_MASS_MSUN:-1e6}"
mkdir -p "${IC_DIR}"

case "${RES}" in
  lowres)
    GRID_IC="${IC_DIR}/agora_lowres_arepo_ic_without_u-with-grid.hdf5"
    OUT_IC="${IC_DIR}/agora_lowres_arepo_ic_with_bh.hdf5"
    EXTRA=(--set-gas-u 0)
    IC_URL="https://www.dropbox.com/scl/fi/m439mkz9vfjvcbflucu3f/agora_lowres_arepo_ic_without_u-with-grid.hdf5?rlkey=c142iij59ruj9nwehj28licsp&dl=1"
    if [ ! -f "${GRID_IC}" ]; then
      echo "get_ics.sh: downloading low-res AGORA IC to ${GRID_IC} ..."
      curl -sL "${IC_URL}" -o "${GRID_IC}"
    fi
    ;;
  medres)
    GRID_IC="${2:-${IC_DIR}/agora_galaxy_med_arepo-with-grid.hdf5}"
    OUT_IC="${IC_DIR}/agora_medres_arepo_ic_with_bh.hdf5"
    EXTRA=(--unit-mass-g 1.989e42 --set-gas-u 0)
    if [ ! -f "${GRID_IC}" ]; then
      echo "get_ics.sh: medium-res IC ${GRID_IC} not found; pass its path as the second argument." >&2
      exit 1
    fi
    ;;
  *)
    echo "usage: $0 [lowres|medres [PATH_TO_IC]]" >&2
    exit 1
    ;;
esac

if [ -f "${OUT_IC}" ]; then
  echo "get_ics.sh: ${OUT_IC} already exists, skipping."
  exit 0
fi

python3 "${HERE}/add_central_bh.py" "${GRID_IC}" "${OUT_IC}" --mass-msun "${BH_MASS_MSUN}" ${EXTRA[@]+"${EXTRA[@]}"}
echo "get_ics.sh: done -> ${OUT_IC}"
