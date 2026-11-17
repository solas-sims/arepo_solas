#!/bin/bash
# Generates the zoom initial conditions for Halo 381 into ./ICs/ by
# running MUSIC2 (https://github.com/cphyc/MUSIC2, Hahn & Abel 2011)
# against MUSIC2_conf/zoom_halo381.conf.
#
# Halo 381 was previously identified by running the L50/h Mpc, N128,
# unigrid parent box (MUSIC2_conf/zoom_parent.conf) to z=0, running FOF
# on it, and tracing the halo's particles back to the parent ICs with
# MUSIC2_scripts/get_zoom_ics_from_parent.py. That analysis produced the
# ref_center/ref_extent Lagrangian patch already baked into
# zoom_halo381.conf -- this script does not redo it, it only runs MUSIC2
# on the already-selected patch. See README.md for the full provenance
# and how to redo the halo selection for a different target.
#
# Set MUSIC2_EXE to point at your MUSIC2 build if it isn't at the
# default sibling-checkout location.

set -euo pipefail

MUSIC2_EXE="${MUSIC2_EXE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../../MUSIC2/build/MUSIC}"

IC_DIR="./ICs"
IC_FILE="${IC_DIR}/L50_N128_z99_Halo381.hdf5"
CONF_FILE="./MUSIC2_conf/zoom_halo381.conf"

mkdir -p "${IC_DIR}"

if [ -f "${IC_FILE}" ]; then
  echo "get_ics.sh: ${IC_FILE} already exists, skipping."
  exit 0
fi

if [ ! -x "${MUSIC2_EXE}" ]; then
  echo "get_ics.sh: MUSIC2 executable not found at '${MUSIC2_EXE}'." >&2
  echo "get_ics.sh: build MUSIC2 (https://github.com/cphyc/MUSIC2) and either" >&2
  echo "get_ics.sh: place the build at that path, or set MUSIC2_EXE to point at it." >&2
  exit 1
fi

echo "get_ics.sh: running MUSIC2 on ${CONF_FILE} ..."
"${MUSIC2_EXE}" "${CONF_FILE}"

if [ ! -f "${IC_FILE}" ]; then
  echo "get_ics.sh: MUSIC2 ran but ${IC_FILE} was not produced -- check zoom_halo381.conf's [output] filename." >&2
  exit 1
fi

if command -v sha256sum >/dev/null 2>&1; then
  sha256sum "${IC_FILE}"
elif command -v shasum >/dev/null 2>&1; then
  shasum -a 256 "${IC_FILE}"
fi

echo "get_ics.sh: done."
