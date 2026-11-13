#!/bin/bash
# Standalone build for the dust rate-law unit test. Deliberately bypasses
# the main Makefile (no MPI/HDF5/GSL needed): the pure rate-law files in
# src/dust/ only depend on dust.h/dust_proto.h and libm.
set -euo pipefail

cd "$(dirname "$0")"

cc -std=c99 -Wall -Wextra -DDUST -o test_dust_rates \
  test_dust_rates.c \
  ../../src/dust/dust_production_sn.c \
  ../../src/dust/dust_destruction_sn.c \
  ../../src/dust/dust_growth.c \
  ../../src/dust/dust_sputtering.c \
  -lm

./test_dust_rates
