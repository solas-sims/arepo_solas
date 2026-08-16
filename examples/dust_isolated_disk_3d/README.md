# dust_isolated_disk_3d (scaffold only, no IC)

This example exists to give the dust module's Phase 1 exit gate ("zero-rate
run reproduces a non-dust run bit-for-bit", "total metal mass conserved,
D/Z-Z trend matches Li et al. 2019") a `Config.sh`/`param.txt` pair to run
against. **It ships with no initial conditions file and cannot be run as-is.**

No example in this repository currently enables `STARS` + `SUPERNOVAE` +
`WINDS` + `METALS` together (the flag combination `DUST` requires) with a
paired IC -- `Config_AGORA.sh` has the right flags but no example directory,
and `examples/galaxy_merger_star_formation_3d/` has an IC but uses the
`EEOS_SF` star-formation model, which never turns on the `STARS`/mechanical-
feedback framework at all. This directory's `Config.sh` mirrors
`Config_AGORA.sh` plus `DUST`; `param.txt` mirrors
`examples/galaxy_merger_star_formation_3d/param.txt`'s general run
parameters plus the extra parameters `AGORA_SF`/`STAR_PARTICLES`/
`STAR_FEEDBACK_ACTIVE`/`SUPERNOVAE`/`METALS`/`USE_GRACKLE` each require.

## What you need to supply before this can run

1. **`InitCondFile`** -- an isolated-disk (or similar idealized ISM) IC.
   Whatever the group's standard isolated-galaxy/ISM-box setup is.
2. **`GrackleDataFile`** -- the Grackle chemistry/cooling data file
   (`CloudyData_UVB=...` etc.), not included in this repo.
3. **`StarTablesFile`** -- the HDF5 stellar yield/feedback table consumed by
   `load_star_tables()` (see `src/stars/star_tables.c`), not included in
   this repo.
4. Adjust `BoxSize`, `ReferenceGasPartMass`, `MeanVolume`, and the
   softening lengths in `param.txt` to match whatever IC you supply --
   the values currently in the file are copied from
   `examples/galaxy_merger_star_formation_3d/param.txt` as placeholders,
   not tuned for this setup.
5. `SN_LeadTime`, `SN_HostShellSweepFrac`, `NumberDensThreshold`,
   `TemperatureThreshold`, `StarFormationEfficiency`, and
   `InitMetallicityinSolar` are all placeholder values -- tune them to
   the group's usual choices for these parameters.

## Running the Phase 1 exit-gate checks once an IC is in place

- **Zero-rate bit-for-bit check**: build once with `DUST` on and once with
  it off (same `Config.sh` otherwise); with the dust rate-law constants in
  `src/dust/dust.h` all effectively zeroed out (e.g.
  `DUST_SN_CONDENSATION_EFFICIENCY 0`), the two runs' gas-cell fields
  (other than `DustMass` itself) should agree to machine precision.
- **Metal mass conservation / D/Z-Z trend**: run with the module's default
  rate-law constants and check `GasMetals` (gas-phase + dust-phase) is
  conserved over the run, and that `DustMass / GasMetals` vs. metallicity
  is qualitatively consistent with Li, Narayanan & Dave (2019)'s published
  relation -- see the "Rate-law coefficient caveat" note in the Phase 1
  implementation plan before treating the latter as a quantitative match.

This example is deliberately **not** added to `test.sh`'s `TESTS` list,
since it cannot run without the inputs above.
