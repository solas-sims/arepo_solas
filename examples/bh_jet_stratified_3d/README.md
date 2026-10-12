# BH jet in a stratified atmosphere (3D)

A smoke test and starting point for `BH_JET_FEEDBACK` with the `BH_JET_LAUNCH_POINT` geometry in the paper's stratified medium (Maragkakis et al. 2026, section 4.5),
scaled down. A black hole at the centre of an isothermal beta-model atmosphere launches bursts that are deposited in gas spheres 5 kpc up and down the z axis. The
atmosphere is held in hydrostatic equilibrium by an external acceleration. With 65,715 cells it runs in 20 seconds on 4 tasks. The uniform-medium version is
`examples/bh_jet_box_3d`, which also explains the jet model in more detail.

It does **not** plug into `test.sh`; run it by hand as described below.

## What is in this directory

| File | Purpose |
|---|---|
| `Config.sh` | Compile-time flags for the jet run |
| `param.txt` | Runtime parameters for the jet run |
| `create.py` | Writes `IC.hdf5`, the beta-model atmosphere with the black hole (needs `numpy` and `h5py`) |
| `check.py` | Post-run checks (energy accounting, momentum, tracer); optional plot |
| `Config_relax.sh`, `param_relax.txt` | Relaxation of the atmosphere without the black hole (optional) |
| `add_bh.py` | Adds the black hole to a relaxed snapshot to make the jet run's ICs (optional) |

## The set-up

* **Atmosphere:** an isothermal beta model, rho = rho_0 [1 + (r/r_c)^2]^(-3 beta / 2) with rho_0 = 1e-25 g/cm^3, beta = 0.38, r_c = 10 kpc and T = 10^7 K (eq. 13). The cells
  have equal mass and are drawn from the profile by inverse-transform sampling inside the sphere that encloses the box, keeping the cells that fall in the box, as the paper
  does for AREPO. The number of cells is set by the first argument to `create.py` (65,715 for 32).
* **Gravity:** no self-gravity. `STATICBETA` (with `EXTERNALGRAVITY`) adds the acceleration that holds the atmosphere in equilibrium, g = -3 beta k T / (mu m_p) **r** / (r^2 + r_c^2),
  centred on the box centre (eq. 14). Without it the gas expands. The constants are in `Config.sh` (`BETA_MODEL_BETA`, `BETA_MODEL_RC`, `BETA_MODEL_TEMP`, `BETA_MODEL_MU`); keep them equal to
  the values in `create.py`. Not for comoving runs. The acceleration also acts on the black hole.
* **Box:** `BoxSize 64` with `LONG_Z=2`, so 64 x 64 x 128 kpc. The hydro box is periodic (the `REFLECTIVE_*` flags are not set, as in the arepobh runs); `GRAVITY_NOT_PERIODIC` and
  `PeriodicBoundariesOn 0` only satisfy a consistency check. The black hole and the jet axis are at the centre and along z.
* **Units:** 1 kpc, 1 km/s, 1e6 Msun, as in the arepobh `param.txt`.
* **Jet:** 1e45 erg/s at 1e4 km/s in bursts of `BhJetBurstMass * v_jet^2` = 2e57 erg (one every ~63 kyr), within 15 degrees of the z axis. Each burst is a pair of virtual particles that
  flies out 5 kpc (0.5 Myr) and is then spread over the gas in a sphere around each end point that holds 8 +- 2 target gas masses. `ENERGY_JET` follows the paper's eq. 12.
  `BH_JET_FIXED_POWER` means the BH's accretion model is irrelevant (`BONDI_ACCRETION` is on only because the build needs one).
* **Run length:** `TimeMax 0.003` is 2.9 Myr.

## 1. Initial conditions

```bash
mkdir -p run && cp Config.sh param.txt run/
python3 create.py run 32
```

The second argument is the number of cells per side in x and y (the box has twice as many in z). `create.py` prints the `ReferenceGasPartMass` that matches (2.57 for 32);
change it in `param.txt` if you change the resolution. `BhJetDesNgb` is in units of that mass.

## 2. Build

From the repository root:

```bash
make -j8 SYSTYPE=MACOSX CONFIG=examples/bh_jet_stratified_3d/run/Config.sh BUILD_DIR=examples/bh_jet_stratified_3d/run/build EXEC=examples/bh_jet_stratified_3d/run/Arepo
```

(`SYSTYPE` as for your machine.) The default flags are `-O0 -g`, about ten times slower. For a quick run add
`OPTIMIZE="-std=c11 -O2 -g -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function"` to the `make` line.

## 3. Run

```bash
cd run
mpirun -np 4 ./Arepo param.txt > log.txt
```

With `-O2` and 4 tasks this takes about 20 seconds and writes five snapshots. Keep the stdout in `log.txt`: `check.py` reads the jet energy accounting from it.
`Config.sh` sets `OVERRIDE_PEANOGRID_WARNING`, which a run on many tasks needs (see `examples/bh_jet_box_3d/README.md`).

## 4. Check

```bash
python3 check.py run
```

It compares the last snapshot with the IC and exits non-zero if a check fails:

1. **Energy.** The change in thermal + kinetic + potential energy of the gas (the potential energy is read from the constants in `Config.sh`) equals the jet energy the code reports as taken
   up by the gas, to 1e-2. In the reference run the difference is 1.2e-3. This is looser than in the uniform example (where it is 1e-13) because the sum also carries the energy errors of the
   gravity coupling and of evaluating the potential at the mesh-generating points. A run without the jet drifts by a similar amount (a few 1e6 over 49 Myr, against the jet's 3.4e9); that
   attribution was not separated more finely. The energy the log reports as given and as taken up by the gas is exactly equal.
2. **Momentum.** The *change* in net momentum since the IC is below 1% of the sum of |p| (1.5e-3 in the reference run). The change is used because a relaxed IC has net momentum already.
3. **Tracer.** Cells carry jet tracer, and they lie around the launch radius and beyond.

It also writes `jet_slice.png`: |v| and the jet tracer in a slab through the BH.

Reference run (4 tasks, `-O2`): 34 bursts delivered by t = 0.003 and 12 still in flight (46 launched), 3.4e9 code units of energy given and taken up, 441 cells with tracer at 3-14 kpc from the BH.
The uniform example delivered 29. The jet energy is produced in whole BH steps, which fall differently in the two runs; this was not examined further.

The energy lines in `log.txt`:

```
BLACKHOLES: Jet energy [code]: produced ..., released ... by 1 of 1 active BHs, held in reservoirs ...
BLACKHOLES: Jet launch points: N packets delivered, ... postponed (empty sphere), ... spheres capped at the launch radius, ... in flight
BLACKHOLES: Energy given by BHs = ..., Energy taken up by gas particles = ...
```

## Equilibrium of the atmosphere

Tested with no jet (`create.py ... nobh`, `Config_relax.sh`, 4 tasks, 196 Myr): the shell densities between 0 and 30 kpc stayed within 3 to 9 per cent of their initial values (ratios to the
initial density at 196 Myr: 1.03, 1.09, 1.04 and 1.00 in the shells 0-5, 5-10, 10-20 and 20-30 kpc) and the mass-weighted rms velocity was about 25 km/s, against a sound speed of about 480 km/s. Most of the
adjustment happens in the first 50 Myr and is the sampling noise of the initial cells. Without the acceleration the central shells fell to 0.33, 0.43, 0.57 and 0.85 of their initial densities in 49 Myr.

## Relaxation

The paper relaxes the sampled atmosphere for 200 Myr with the potential on before launching the jet, which removes the sampling noise.

1. `python3 create.py relax 32 nobh`, copy `Config_relax.sh` and `param_relax.txt` into `relax/` as `Config.sh` and `param.txt`, build and run as above. It writes snapshots every 49 Myr;
   `snap_004.hdf5` is at t = 0.2 (196 Myr). No black hole and no jet. This takes about 9 minutes at `-O2` on 4 tasks.
2. `python3 add_bh.py relax/output/snap_004.hdf5 run/IC.hdf5` writes the jet run's ICs: the relaxed gas as it is (including its leftover motions) plus the black hole at the box centre,
   with the header time reset to 0. The BH mass is an optional third argument (default 1, which is 1e6 Msun). `LONG_Z` in the script must match `Config.sh`.
3. Run the jet as above.

This chain was run at the 32-cell resolution (a 196 Myr relaxation to 65,708 cells, then the jet). The relaxed gas keeps random motions of about 25 km/s rms.

## Switching the jet off

For the paper's remnant runs. Add `BH_JET_STOP_TIME` to `Config.sh` (a commented line is there) and `BhJetStopTime` to `param.txt` (code time; 0.045 is 44 Myr). From that time on no jet energy is produced and no
burst is launched, and the unreleased energy is dropped. Bursts already in flight are still delivered. In a test on the relaxed atmosphere with `BhJetStopTime 0.0015` in a run to 0.003, bursts were launched
at t = 0.00075 and 0.001125 (17 in all), the last was delivered at 0.0019, and nothing was produced after 0.0015; the energy given and taken up by the gas were equal (1.7e9, 17 bursts).

The paper's remnant runs (its section 4.6) use a 900 x 900 x 1200 kpc box with 1.8e7 elements of 10^6 Msun, a 44 Myr jet with `m_j = 2e4 Msun` (`BhJetBurstMass 2e4`), and 196 Myr in all
(`TimeMax 0.2`). That would need `BoxSize 900` and `LONG_Z` = 4/3, a non-integer value that has not been tried, and it has not been run.

## Relation to the paper

The paper's stratified runs (its section 4.5) differ from this example as follows.

| | Paper | This example |
|---|---|---|
| Atmosphere | rho_0 = 1e-25 g/cm^3, beta = 0.38, r_c = 10 kpc, T = 1e7 K | same |
| Box | 300 x 300 x 600 kpc, periodic hydro | 64 x 64 x 128 kpc |
| Elements | N = 2.9e7 of 1e5 Msun, sampled in a sphere, relaxed 200 Myr | 65,715 of 2.6e6 Msun, sampled in a sphere, not relaxed unless you follow the relaxation steps |
| Jet power, speed | 5e45 erg/s, 4e4 km/s | 1e45, 1e4 |
| Burst mass | `m_j = 2e4 Msun` (halved from the uniform runs) | `BhJetBurstMass 1e6` Msun |
| Duration | 44 Myr | 2.9 Myr |
| Launch radius, cone | 10 kpc, 15 degrees | 5 kpc, 15 degrees |
| Deposit region | the 10 nearest gas elements | a sphere holding 8 +- 2 target gas masses (the paper's section 4.4 scheme) |

To get closer, set `BoxSize 300`, `TimeMax 0.045`, `BhJetPower 5e45`, `BhJetVelocity 4e4`, `BhJetBurstMass 2e4`, `BhJetLaunchRadius 10`, `BhJetDesNgb 10` and `ReferenceGasPartMass 0.1`, use a much larger grid in
`create.py`, and relax first. That is a cluster-sized run. In the paper's results for this set-up (its Fig. 12) the lobes are longer and thinner than in the uniform medium, and the jets heat the gas in the inner regions and displace it to larger radii.

## Not tested here

Comoving integration (`STATICBETA` is non-comoving only), the remnant box, a BH that moves between tasks while bursts are in flight, and anything at paper scale.
