#!/bin/bash        	# this line only there to enable syntax highlighting in this file

#########################################################
#  Enable/Disable compile-time options as needed        #
#  examples/galaxy_merger_star_formation_3d/Config.sh   #
#########################################################

#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT                # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE              # Use maximum face angle as roundness criterion in mesh regularization

#--------------------------------------- Refinement and derefinement
REFINEMENT_SPLIT_CELLS 			 # Refinement
REFINEMENT_MERGE_CELLS                   # Derefinement
#REFINEMENT_VOLUME_LIMIT    		 # Limit the volume of cells and the maximum volume difference between neighboring cels
#NODEREFINE_BACKGROUND_GRID 		 # Do not de-refine low-res gas cells in zoom simulations

#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS 			 # non-local timestep criterion (take 'signal speed' into account)

#------------------------------------------------ TreePM Options
PMGRID=256                                       # Enables particle mesh; numberof cells used for grid in each dimension
RCUT=5.0                                         # This can be used to override the maximum radius in which the short-range tree-force is evaluated (in case the TreePM algorithm is used). The default value is 4.5, given in mesh-cells.

#--------------------------------------- Gravity treatment
SELFGRAVITY                    # gravitational intraction between simulation particles/cells
HIERARCHICAL_GRAVITY           # use hierarchical splitting of the time integration of the gravity
CELL_CENTER_GRAVITY            # uses geometric centers to calculate gravity of cells, only possible with HIERARCHICAL_GRAVITY
ALLOW_DIRECT_SUMMATION         # Performed direct summation instead of tree-based gravity if number of active particles < DIRECT_SUMMATION_THRESHOLD (= 3000 unless specified differently here)
DIRECT_SUMMATION_THRESHOLD=500 # Overrides maximum number of active particles for which direct summation is performed instead of tree based calculation

#--------------------------------------- Gravity softening
NSOFTTYPES=4                    # Number of different softening values to which particle types can be mapped.
MULTIPLE_NODE_SOFTENING         # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
#INDIVIDUAL_GRAVITY_SOFTENING=32 # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
ADAPTIVE_HYDRO_SOFTENING        # Adaptive softening of gas cells depending on their size

#--------------------------------------- Single/Double Precision
DOUBLEPRECISION=1        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
NGB_TREE_DOUBLEPRECISION # if this is enabled, double precision is used for the neighbor node extension

#-------------------------------------------- Things for special behaviour
PROCESS_TIMES_OF_OUTPUTLIST # goes through times of output list prior to starting the simulaiton to ensure that outputs are written as close to the desired time as possible (as opposed to at next possible time if this flag is not active)
OVERRIDE_PEANOGRID_WARNING  # don't stop if peanogrid is not fine enough

#--------------------------------------- Output/Input options
HAVE_HDF5 # needed when HDF5 I/O support is desired (recommended)
INPUT_IN_DOUBLEPRECISION
GENERATE_GAS_IN_ICS

#------------------------------------------------ On the fly FOF groupfinder
FOF                                              # enable FoF output
FOF_PRIMARY_LINK_TYPES=2                         # 2^type for the primary dark matter type
FOF_SECONDARY_LINK_TYPES=1+16+32                 # 2^type for the types linked to nearest primaries

#------------------------------------------------ Subfind
SUBFIND                                          # enables substructure finder
#SAVE_HSML_IN_SNAPSHOT                            # stores hsml, density, and velocity dispersion values in the snapshot files
#SUBFIND_CALC_MORE                                # calculates also the velocity dispersion in the local density estimate (this is automatically enabled by several other options, e.g. SAVE_HSML_IN_SNAPSHOT)
#SUBFIND_EXTENDED_PROPERTIES                      # adds calculation of further quantities related to angular momentum in different components

#--------------------------------------- Testing and Debugging options
DEBUG # enables core-dumps

#--------------------------------------- non-standard phyiscs
ENFORCE_JEANS_STABILITY_OF_CELLS # this imposes an adaptive floor for the temperature
COOLING                          # Simple primordial cooling
USE_GRACKLE			 # Use Grackle library for cooling
GRACKLE_CHEMISTRY=0
#NOUVBACKGROUND			 # Switches off UV background
USE_SFR                          # Activate star formation
#EEOS_SF			 # AGORA star formation scheme
JEANS_SF
#JEANS_MASS_BASED
USE_CELIB
STARS
SUPERNOVAE
#PASSIVE_SCALARS=1
#METALS

#REFINEMENT_AROUND_BH
#BLACKHOLES
#BONDI_ACCRETION
