#!/bin/bash
# Downloads the two external data tables this example needs but doesn't
# ship: a Grackle cooling table (Cloudy metals 2008, 3D) and the stellar feedback lookup table
# read by src/stars/star_tables.c. See README.md for what these are and
# the schema the star tables file must follow.

set -euo pipefail

GRACKLE_DIR="./grackle_data"
GRACKLE_FILE="${GRACKLE_DIR}/cloudy_metals_2008_3D.h5"
GRACKLE_URL="https://www.dropbox.com/scl/fi/ggiiqlgucb4dkm4louubj/cloudy_metals_2008_3D.h5?rlkey=x9b9amouxqjfu1m023e8gyv31&dl=1"

STAR_TABLES_DIR="./star_tables"
STAR_TABLES_FILE="${STAR_TABLES_DIR}/star_feedback_tables.hdf5"
STAR_TABLES_URL="https://www.dropbox.com/scl/fi/a085tatarnwrx3we553eu/star_feedback_tables.hdf5?rlkey=t07hdi3t5k18868kojd6nsk1y&dl=1"

mkdir -p "${GRACKLE_DIR}" "${STAR_TABLES_DIR}"

if [ -f "${GRACKLE_FILE}" ]; then
  echo "get_tables.sh: ${GRACKLE_FILE} already exists, skipping download."
else
  echo "get_tables.sh: downloading Grackle cooling table to ${GRACKLE_FILE} ..."
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
