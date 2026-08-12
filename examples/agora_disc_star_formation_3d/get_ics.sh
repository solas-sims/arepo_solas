#!/bin/bash
# Downloads the AGORA isolated-disc initial conditions (low-res, Arepo HDF5
# format, already passed through the ADDBACKGROUNDGRID step) into ./ICs/.
#
# Source: Dropbox share provided by the SOLAS group.
# SHA256 of the file as originally verified when this script was written:
#   4f5d58aa9d1d7e91862d0a7f1d359e0c75e8e70f681ddbc425a6c1f15cc063f
# (informational only -- not enforced, since the upstream file may be
# replaced with an updated version).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IC_DIR="${SCRIPT_DIR}/ICs"
IC_FILE="${IC_DIR}/agora_lowres_arepo_ic_without_u-with-grid.hdf5"
IC_URL="https://www.dropbox.com/scl/fi/m439mkz9vfjvcbflucu3f/agora_lowres_arepo_ic_without_u-with-grid.hdf5?rlkey=c142iij59ruj9nwehj28licsp&dl=1"

mkdir -p "${IC_DIR}"

if [ -f "${IC_FILE}" ]; then
  echo "get_ics.sh: ${IC_FILE} already exists, skipping download."
  exit 0
fi

echo "get_ics.sh: downloading AGORA disc ICs to ${IC_FILE} ..."
curl -sL "${IC_URL}" -o "${IC_FILE}"

if command -v sha256sum >/dev/null 2>&1; then
  sha256sum "${IC_FILE}"
elif command -v shasum >/dev/null 2>&1; then
  shasum -a 256 "${IC_FILE}"
fi

echo "get_ics.sh: done."
