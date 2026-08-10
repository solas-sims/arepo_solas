#--------------------------------------- SOLAS additions

#--------------------------------------- Metal parameters
METALS                 # Advect all metals, ie metal mass fraction, as a PASSIVE_SCALARS

#--------------------------------------- Cooling parameters
USE_GRACKLE
GRACKLE_CHEMISTRY=3    

#--------------------------------------- Star Formation options
AGORA_SF               # Agora based SF

#--------------------------------------- Star options
STARS                  # General stars framework flag

STAR_PARTICLES=1       # Star particles model flag: set to 0, 1 for massive star particles, set to 2 for resolved individual stars

#STAR_FEEDBACK          # Include full star feedback (winds + full radiation + supernovae)
WINDS                  # Only winds
#RADIATION              # Full radiation
SUPERNOVAE             # Only supernovae

#STAR_HOST_REFINEMENT


#--------------------------------------- Arepo public

#--------------------------------------- Mesh motion and regularization; default: moving mesh
REGULARIZE_MESH_CM_DRIFT      # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE    # Use maximum face angle as roundness criterion in mesh regularization

#--------------------------------------- Refinement and derefinement; default: no refinement/derefinement; criterion: target mass
REFINEMENT_SPLIT_CELLS        # Refinement
REFINEMENT_MERGE_CELLS        # Derefinement
REFINEMENT_VOLUME_LIMIT       # Limit the volume of cells and the maximum volume difference between neighboring cels
NODEREFINE_BACKGROUND_GRID    # Do not de-refine low-res gas cells in zoom simulations

#--------------------------------------- non-standard phyiscs
COOLING                       # Simple primordial cooling
USE_SFR                       # Star formation model, turning dense gas into collisionless partices

#--------------------------------------- Gravity treatment; default: no gravity
SELFGRAVITY                   # gravitational intraction between simulation particles/cells
HIERARCHICAL_GRAVITY          # use hierarchical splitting of the time integration of the gravity
CELL_CENTER_GRAVITY           # uses geometric centers to calculate gravity of cells, only possible with HIERARCHICAL_GRAVITY
GRAVITY_NOT_PERIODIC          # gravity is not treated periodically
ALLOW_DIRECT_SUMMATION        # Performed direct summation instead of tree-based gravity if number of active particles < DIRECT_SUMMATION_THRESHOLD (= 3000 unless specified differently here)
DIRECT_SUMMATION_THRESHOLD=500  # Overrides maximum number of active particles for which direct summation is performed instead of tree based calculation

#--------------------------------------- Gravity softening
NSOFTTYPES=4                  # Number of different softening values to which particle types can be mapped.
MULTIPLE_NODE_SOFTENING       # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
ADAPTIVE_HYDRO_SOFTENING      # Adaptive softening of gas cells depending on their size

#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS          # non-local timestep criterion (take 'signal speed' into account)

#--------------------------------------- Single/Double Precision
DOUBLEPRECISION=1             # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
NGB_TREE_DOUBLEPRECISION      # if this is enabled, double precision is used for the neighbor node extension

#--------------------------------------- output options
PROCESS_TIMES_OF_OUTPUTLIST   # goes through times of output list prior to starting the simulaiton to ensure that outputs are written as close to the desired time as possible (as opposed to at next possible time if this flag is not active)
HAVE_HDF5                     # needed when HDF5 I/O support is desired (recommended)

#--------------------------------------- Testing and Debugging options
DEBUG                         # enables core-dumps

OVERRIDE_PEANOGRID_WARNING  # don't stop if peanogrid is not fine enough