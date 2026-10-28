#!/bin/bash
# Prepares the AGORA isolated-disc initial conditions *with a central black hole*
# (1e6 Msun PartType5 at the disc centre) in ./ICs/.
#
#   ./get_ics.sh lowres                 downloads the low-res IC and adds the BH
#   ./get_ics.sh medres [PATH_TO_IC]    downloads the medium-res IC (or uses PATH_TO_IC) and adds the BH
#
# Needs python3 with h5py and numpy. BH mass: BH_MASS_MSUN (default 1e6).
#
# The medium-res IC is the AGORA MED set already passed through ADDBACKGROUNDGRID.
# Both ICs already have the right UnitMass_in_g (1.989e42) and a gas InternalEnergy of ~101 (code units),
# so the only change made is adding the BH at BoxSize/2, which is the disc centre in both files.

set -euo pipefail

RES="${1:-lowres}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IC_DIR="./ICs"
BH_MASS_MSUN="${BH_MASS_MSUN:-1e6}"
mkdir -p "${IC_DIR}"

case "${RES}" in
  lowres)
    GRID_IC="${IC_DIR}/agora_lowres_arepo_ic-with-grid.hdf5"
    OUT_IC="${IC_DIR}/agora_lowres_arepo_ic_with_bh.hdf5"
    EXTRA=()
    IC_URL="https://www.dropbox.com/scl/fi/86uaihmoairvnpuexe5ad/agora_lowres_arepo_ic-with-grid.hdf5?rlkey=j9o95uah84tgeplb9k8ladxu2&dl=1"
    if [ ! -f "${GRID_IC}" ]; then
      echo "get_ics.sh: downloading low-res AGORA IC to ${GRID_IC} ..."
      curl -sL "${IC_URL}" -o "${GRID_IC}"
    fi
    ;;
  medres)
    GRID_IC="${2:-${IC_DIR}/agora_medres_arepo_ic-with-grid.hdf5}"
    OUT_IC="${IC_DIR}/agora_medres_arepo_ic_with_bh.hdf5"
    EXTRA=()
    IC_URL="https://www.dropbox.com/scl/fi/f5h135z2bkpcu4axp3y4c/agora_medres_arepo_ic-with-grid.hdf5?rlkey=yj1esq1ug4v5oulfka6e9re1h&dl=1"
    if [ ! -f "${GRID_IC}" ]; then
      echo "get_ics.sh: downloading medium-res AGORA IC to ${GRID_IC} ..."
      curl -sL "${IC_URL}" -o "${GRID_IC}"
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
