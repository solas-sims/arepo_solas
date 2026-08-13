#!/bin/bash
# Downloads the L10N128 cosmological initial conditions (dark-matter-only,
# 128^3 particles in a 10 Mpc/h box; gas is split off at startup via
# GENERATE_GAS_IN_ICS) into ./ICs/.
#
# Source: Dropbox share provided by the SOLAS group.

set -euo pipefail

IC_DIR="./ICs"
IC_FILE="${IC_DIR}/cosmo_L10N128.hdf5"
IC_URL="https://www.dropbox.com/scl/fi/rqi9a78390t2cb2sliugt/cosmo_L10N128.hdf5?rlkey=7l98ph0tj3fm3fm0mr6oncofy&dl=1"

mkdir -p "${IC_DIR}"

if [ -f "${IC_FILE}" ]; then
  echo "get_ics.sh: ${IC_FILE} already exists, skipping download."
  exit 0
fi

echo "get_ics.sh: downloading L10N128 cosmo ICs to ${IC_FILE} ..."
curl -sL "${IC_URL}" -o "${IC_FILE}"

if command -v sha256sum >/dev/null 2>&1; then
  sha256sum "${IC_FILE}"
elif command -v shasum >/dev/null 2>&1; then
  shasum -a 256 "${IC_FILE}"
fi

echo "get_ics.sh: done."
