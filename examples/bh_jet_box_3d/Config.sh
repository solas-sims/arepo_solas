#--------------------------------------- BH jet in a hot-gas box (Maragkakis et al. 2026 set-up, scaled down)
#
# Periodic hydro box, stretched in z (BoxSize x BoxSize x 2*BoxSize), gravity off.
# The REFLECTIVE_* flags are deliberately NOT set: the hydro box is periodic, as in the arepobh runs.
# GRAVITY_NOT_PERIODIC is set (so param.txt has PeriodicBoundariesOn 0); with no SELFGRAVITY it only
# satisfies the consistency check in check_parameters().

#--------------------------------------- Box
LONG_Z=2                       # stretch the z extent of the box by this (integer) factor
GRAVITY_NOT_PERIODIC           # gravity is not treated periodically

#--------------------------------------- Black hole and jet
BLACKHOLES
BONDI_ACCRETION                # required by BH_FEEDBACK_ACTIVE; irrelevant here (BH_JET_FIXED_POWER ignores the accretion rate)
BH_JET_FEEDBACK                # jet injection (src/blackholes/bh_jet_*.c)
ENERGY_JET                     # inject the jet energy as kinetic energy, conserving total energy exactly
BH_JET_FIXED_POWER             # fixed power BhJetPower [erg/s] instead of Epsilon_f Epsilon_r Mdot c^2; creates the jet mass
BH_JET_LAUNCH_POINT            # arepobh geometry: bursts are deposited in gas spheres at BhJetLaunchRadius along +-axis
BURST_MODE                     # release the energy in bursts of BhJetBurstMass * v_jet^2
JET_TRACER                     # passive tracer = 1 in cells that receive jet material; such cells are not de-refined
#BH_JET_PRESSURE_TRIGGER       # optional: only launch if the jet ram pressure exceeds BhJetPressureFactor * ambient pressure
#BH_JET_SPIN                   # optional: jet axis follows the accreted angular momentum instead of BhJetAxisX/Y/Z

#--------------------------------------- Mesh
REGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization
REFINEMENT_SPLIT_CELLS                   # Refinement
REFINEMENT_MERGE_CELLS                   # Derefinement
TREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)

#--------------------------------------- Precision, I/O
DOUBLEPRECISION=1              # full double precision
INPUT_IN_DOUBLEPRECISION       # initial conditions are in double precision
OUTPUT_IN_DOUBLEPRECISION      # snapshot files will be written in double precision
HAVE_HDF5                      # HDF5 I/O
DEBUG                          # enables core-dumps
OVERRIDE_PEANOGRID_WARNING     # don't stop if the top-tree cannot be split further (see README: needed with many tasks)
