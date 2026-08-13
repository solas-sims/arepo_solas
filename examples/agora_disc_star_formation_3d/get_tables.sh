#!/bin/bash
# Downloads the two external data tables this example needs but doesn't
# ship: a Grackle cooling/UVB table and the stellar feedback lookup table
# read by src/stars/star_tables.c. See README.md for what these are and
# the schema the star tables file must follow.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

GRACKLE_DIR="${SCRIPT_DIR}/grackle_data"
GRACKLE_FILE="${GRACKLE_DIR}/CloudyData_UVB=FG2011_shielded.h5"
GRACKLE_URL="https://raw.githubusercontent.com/grackle-project/grackle_data_files/928696482fbe15d9bac4382de6134d95568f099c/input/CloudyData_UVB%3DFG2011_shielded.h5"

STAR_TABLES_DIR="${SCRIPT_DIR}/star_tables"
STAR_TABLES_FILE="${STAR_TABLES_DIR}/star_feedback_tables.hdf5"
STAR_TABLES_URL="https://www.dropbox.com/scl/fi/a085tatarnwrx3we553eu/star_feedback_tables.hdf5?rlkey=t07hdi3t5k18868kojd6nsk1y&dl=1"

mkdir -p "${GRACKLE_DIR}" "${STAR_TABLES_DIR}"

if [ -f "${GRACKLE_FILE}" ]; then
  echo "get_tables.sh: ${GRACKLE_FILE} already exists, skipping download."
else
  echo "get_tables.sh: downloading Grackle UVB table to ${GRACKLE_FILE} ..."
  curl -sL "${GRACKLE_URL}" -o "${GRACKLE_FILE}"
fi

if [ -f "${STAR_TABLES_FILE}" ]; then
  echo "get_tables.sh: ${STAR_TABLES_FILE} already exists, skipping download."
else
  echo "get_tables.sh: downloading stellar feedback tables to ${STAR_TABLES_FILE} (~180 MB) ..."
  curl -sL "${STAR_TABLES_URL}" -o "${STAR_TABLES_FILE}"
fi

if command -v sha256sum >/dev/null 2>&1; then
  sha256sum "${GRACKLE_FILE}" "${STAR_TABLES_FILE}"
elif command -v shasum >/dev/null 2>&1; then
  shasum -a 256 "${GRACKLE_FILE}" "${STAR_TABLES_FILE}"
fi

echo "get_tables.sh: done."
