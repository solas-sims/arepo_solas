# AREPO Makefile
EXEC = Arepo
LIBRARY = arepo
CONFIG = Config.sh
BUILD_DIR = build
SRC_DIR = src

SYSTYPE ?= "LINUX"

MAKEFILES = Makefile config-makefile
ifeq ($(wildcard Makefile.systype), Makefile.systype)
MAKEFILES += Makefile.systype
endif

include config-makefile
-include Makefile.systype

$(info Build configuration:)
$(info SYSTYPE: $(SYSTYPE))
$(info CONFIG: $(CONFIG))
$(info EXEC: $(EXEC))
$(info )

PYTHON = python3
PERL = /usr/bin/perl
RESULT := $(shell CONFIG=$(CONFIG) PERL=$(PERL) BUILD_DIR=$(BUILD_DIR) make -f config-makefile)
CONFIGVARS := $(shell cat $(BUILD_DIR)/arepoconfig.h)
RESULT := $(shell SRC_DIR=$(SRC_DIR) BUILD_DIR=$(BUILD_DIR) ./git_version.sh)

# Default
MPICH_INCL =
MPICH_LIB  = -lmpich
GMP_LIB    = -lgmp
GSL_LIB    = -lgsl -lgslcblas
MATH_LIB   = -lm -lstdc++
HWLOC_LIB  = -lhwloc

################################
#determine the needed libraries#
################################

# we only need fftw if PMGRID is turned on, make sure variable is empty otherwise
FFTW_LIB =
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
FFTW_LIB = $(FFTW_LIBS) -lfftw3
else
FFTW_LIB = $(FFTW_LIBS) -lfftw3f
endif
endif

ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(CONFIGVARS)))
HDF5_INCL = 
HDF5_LIB = 
endif

ifneq (IMPOSE_PINNING,$(findstring IMPOSE_PINNING,$(CONFIGVARS)))
HWLOC_INCL =
HWLOC_LIB =
endif

ifeq (USE_GRACKLE,$(findstring USE_GRACKLE,$(CONFIGVARS)))
OPTIONS += -DCONFIG_BFLOAT_8
GRACKLE_INCL = -I$(HOME)/Codes/grackle/include
GRACKLE_LIB = -L$(HOME)/Codes/grackle/lib -lgrackle -Wl,-rpath,$(HOME)/Codes/grackle/lib
ifeq ($(SYSTYPE),"MACOSX")
LDFLAGS += -L$(shell BREW --prefix gcc)/lib/gcc/current -lgfortran -lquadmath
else
LDFLAGS += -lgfortran -lquadmath
endif
else
GRACKLE_INCL =
GRACKLE_LIB =
LDFLAGS +=
endif

#Mac OS using MacPorts modules for openmpi, fftw, gsl, hdf5 and hwloc
ifeq ($(filter Darwin,$(SYSTYPE)),Darwin)
# compiler and its optimization options
CC        = mpicc   # sets the C-compiler
OPTIMIZE  = -std=c11 -ggdb -g -O0 -fno-omit-frame-pointer -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function

MPICH_LIB = -lmpi
GSL_INCL  = -I/opt/local/include
GSL_LIB   = -L/opt/local/lib -lgsl -lgslcblas
HWLOC_LIB = -L/opt/local/lib -lhwloc

# libraries that are included on demand, depending on Config.sh options
FFTW_INCL = -I/opt/local/include -I/usr/local/include
FFTW_LIBS = -L/opt/local/lib -I/usr/local/lib
HDF5_INCL = -I/opt/local/include 
HDF5_LIB  = -L/opt/local/lib  
HWLOC_INCL= -I/opt/local/include
endif
# end of Darwin

#Mac OS using MacPorts modules for openmpi, fftw, gsl, hdf5 and hwloc
ifeq ($(filter MACOSX,$(SYSTYPE)),MACOSX)
BREW := /opt/homebrew/bin/brew
$(info BREW: $(BREW))

# compiler and its optimization options
CC        = mpicc   # sets the C-compiler
OPTIMIZE  = -std=c11 -ggdb -g -O0 -fno-omit-frame-pointer -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function

MPICH_LIB = #-lmpi
GSL_INCL  = -I$(shell $(BREW) --prefix gsl)/include
GSL_LIB   = -L$(shell $(BREW) --prefix gsl)/lib -lgsl -lgslcblas
GMP_INCL  = -I$(shell $(BREW) --prefix gmp)/include
GMP_LIB   = -L$(shell $(BREW) --prefix gmp)/lib -lgmp

# libraries that are included on demand, depending on Config.sh options
FFTW_INCL = -I$(shell $(BREW) --prefix fftw)/include -I/usr/local/include
FFTW_LIBS = -L$(shell $(BREW) --prefix fftw)/lib -I/usr/local/lib
HDF5_INCL = -I$(shell $(BREW) --prefix hdf5)/include 
HDF5_LIB  = -L$(shell $(BREW) --prefix hdf5)/lib 
HWLOC_INCL= -I$(shell $(BREW) --prefix hwloc)/include
HWLOC_LIB = -L$(shell $(BREW) --prefix hwloc)/lib -lhwloc
endif
# end of Darwin

#Linux
ifeq ($(filter LINUX,$(SYSTYPE)),"LINUX")
# compiler and its optimization options
CC        = mpicc
OPTIMIZE  = -std=c11 -ggdb -g -O0 -fno-omit-frame-pointer -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function

MPICH_INCL= -I/usr/lib/x86_64-linux-gnu/openmpi/include/
MPICH_LIB = -L/usr/lib/x86_64-linux-gnu/openmpi/lib/ -lmpi
GSL_INCL  =
GSL_LIB   = -lgsl -lgslcblas
HWLOC_LIB = -lhwloc

# libraries that are included on demand, depending on Config.sh options
FFTW_INCL =
FFTW_LIBS =
HDF5_INCL = -I/usr/include/hdf5/serial/ 
HDF5_LIB  = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/
HWLOC_INCL=
endif
# end of Linux

#Ngarrgu Tindebeek
ifeq ($(filter NT,$(SYSTYPE)),"NT")
# compiler and its optimization options
CC        =  mpicc
OPTIMIZE  =  -std=c11 -ggdb -O3 -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function

MPICH_INCL=
MPICH_LIB = -lmpi
GSL_INCL  = -I$(EBROOTGSL)/include
GSL_LIB   = -L$(EBROOTGSL)/lib -lgsl -lgslcblas
HWLOC_LIB =

# libraries that are included on demand, depending on Config.sh options
FFTW_INCL = -I$(EBROOTFFTW)/include
FFTW_LIBS = -L$(EBROOTFFTW)/lib
HDF5_INCL = -I$(EBROOTHDF5)/include -DH5_USE_16_API
HDF5_LIB  = -L$(EBROOTHDF5)/lib -lhdf5 -lz
HWLOC_INCL=
endif
# end of NT

# Pawsey Setonix (Cray: PrgEnv-cray + cray-mpich; load gsl, fftw, hdf5 modules)
ifneq ($(filter Setonix "Setonix",$(SYSTYPE)),)
# the Cray compiler wrapper 'cc' provides MPI includes and libraries itself
CC        = cc
OPTIMIZE  = -std=c11 -ggdb -g -O0 -fno-omit-frame-pointer -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function
#OPTIMIZE = -std=c11 -g -O2 -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function  # production

MPICH_INCL=
MPICH_LIB = #-lmpi provided by the cc wrapper

GSL_INCL  = -I$(PAWSEY_GSL_HOME)/include
GSL_LIB   = -L$(PAWSEY_GSL_HOME)/lib -lgsl -lgslcblas

# NOTE: keep gmp in its own prefix. A shared prefix (e.g. ~/software) whose
# include/ also contains another MPI's mpi.h will shadow cray-mpich's mpi.h
# and break the link with undefined ompi_* symbols.
GMP_INCL  = -I$(HOME)/software/gmp-6.3.0/include
GMP_LIB   = -L$(HOME)/software/gmp-6.3.0/lib64 -lgmp

# libraries that are included on demand, depending on Config.sh options
FFTW_INCL = -I$(PAWSEY_FFTW_HOME)/include
FFTW_LIBS = -L$(PAWSEY_FFTW_HOME)/lib
HDF5_INCL = -I$(PAWSEY_HDF5_HOME)/include -DH5_USE_16_API
HDF5_LIB  = -L$(PAWSEY_HDF5_HOME)/lib -lhdf5 -lz
HWLOC_INCL=
HWLOC_LIB =
endif
# end of Setonix

ifndef LINKER
LINKER = $(CC)
endif

##########################################
#determine the needed object/header files#
##########################################

SUBDIRS = . \
          debug_md5 \
          domain \
          gitversion \
          gravity \
          gravity/pm \
          hydro \
          init \
          io \
          main \
          mesh \
          mesh/voronoi \
          mpi_utils \
          ngbtree \
          pm \
          time_integration \
          utils \

OBJS = debug_md5/calc_checksum.o \
       debug_md5/Md5.o \
       domain/domain.o \
       domain/domain_balance.o \
       domain/domain_box.o \
       domain/domain_counttogo.o \
       domain/domain_DC_update.o \
       domain/domain_exchange.o \
       domain/domain_rearrange.o \
       domain/domain_sort_kernels.o \
       domain/domain_toplevel.o \
       domain/domain_vars.o \
       domain/peano.o \
       gravity/accel.o \
       gravity/forcetree.o \
       gravity/forcetree_ewald.o  \
       gravity/forcetree_optimizebalance.o \
       gravity/forcetree_walk.o \
       gravity/grav_external.o \
       gravity/grav_softening.o \
       gravity/gravdirect.o \
       gravity/gravtree.o \
       gravity/gravtree_forcetest.o \
       gravity/longrange.o \
       gravity/pm/pm_periodic2d.o \
       gravity/pm/pm_periodic.o \
       gravity/pm/pm_mpi_fft.o \
       gravity/pm/pm_nonperiodic.o \
       hydro/finite_volume_solver.o \
       hydro/gradients.o \
       hydro/riemann.o \
       hydro/riemann_hllc.o \
       hydro/riemann_hlld.o \
       hydro/scalars.o \
       hydro/update_primitive_variables.o \
       init/begrun.o \
       init/density.o \
       init/init.o \
       io/global.o \
       io/hdf5_util.o \
       io/io.o \
       io/io_fields.o \
       io/logs.o \
       io/parameters.o \
       io/read_ic.o \
       io/restart.o \
       main/allvars.o \
       main/main.o \
       main/run.o \
       mesh/criterion_derefinement.o \
       mesh/criterion_refinement.o \
       mesh/refinement.o \
       mesh/set_vertex_velocities.o \
       mesh/voronoi/voronoi.o \
       mesh/voronoi/voronoi_1d.o \
       mesh/voronoi/voronoi_1d_spherical.o \
       mesh/voronoi/voronoi_3d.o \
       mesh/voronoi/voronoi_check.o \
       mesh/voronoi/voronoi_derefinement.o \
       mesh/voronoi/voronoi_dynamic_update.o \
       mesh/voronoi/voronoi_exchange.o \
       mesh/voronoi/voronoi_ghost_search.o \
       mesh/voronoi/voronoi_gradients_lsf.o \
       mesh/voronoi/voronoi_gradients_onedims.o \
       mesh/voronoi/voronoi_refinement.o \
       mesh/voronoi/voronoi_utils.o \
       mpi_utils/checksummed_sendrecv.o \
       mpi_utils/hypercube_allgatherv.o \
       mpi_utils/mpi_util.o \
       mpi_utils/myalltoall.o \
       mpi_utils/sizelimited_sendrecv.o \
       mpi_utils/pinning.o \
       ngbtree/ngbtree.o \
       ngbtree/ngbtree_search.o \
       ngbtree/ngbtree_walk.o \
       time_integration/darkenergy.o \
       time_integration/do_gravity_hydro.o \
       time_integration/driftfac.o \
       time_integration/predict.o \
       time_integration/timestep.o \
       time_integration/timestep_treebased.o \
       utils/allocate.o \
       utils/debug.o \
       utils/mpz_extension.o \
       utils/mymalloc.o \
       utils/parallel_sort.o \
       utils/system.o \

INCL += debug_md5/Md5.h \
        domain/bsd_tree.h \
        domain/domain.h \
        gitversion/version.h\
        gravity/forcetree.h \
        main/allvars.h \
        main/proto.h \
        mesh/mesh.h \
        mesh/voronoi/voronoi.h \
        time_integration/timestep.h \
        utils/dtypes.h \
        utils/generic_comm_helpers2.h \
        utils/timer.h

ifeq (TWODIMS,$(findstring TWODIMS,$(CONFIGVARS)))
OBJS += mesh/voronoi/voronoi_2d.o
endif

ifeq (MYIBARRIER,$(findstring MYIBARRIER,$(CONFIGVARS)))
OBJS += mpi_utils/myIBarrier.o
INCL += mpi_utils/myIBarrier.h
endif

ifeq (MHD,$(findstring MHD,$(CONFIGVARS)))
OBJS += hydro/mhd.o
endif

ifeq (ADDBACKGROUNDGRID,$(findstring ADDBACKGROUNDGRID,$(CONFIGVARS)))
OBJS += add_backgroundgrid/add_bggrid.o \
        add_backgroundgrid/calc_weights.o \
        add_backgroundgrid/distribute.o
INCL += add_backgroundgrid/add_bggrid.h
SUBDIRS += add_backgroundgrid
endif

#COOLING
ifeq (COOLING,$(findstring COOLING,$(CONFIGVARS)))
OBJS += cooling/cooling.o
INCL += cooling/cooling_vars.h \
        cooling/cooling_proto.h
SUBDIRS += cooling
endif

ifeq (USE_GRACKLE,$(findstring USE_GRACKLE,$(CONFIGVARS)))
OBJS += cooling/grackle.o
endif

#SFR
ifneq (,$(filter USE_SFR,$(CONFIGVARS)))
OBJS += star_formation/starformation.o 
SUBDIRS += star_formation 

ifneq (,$(filter EEOS_SF AGORA_SF JEANS_SF,$(CONFIGVARS)))
ifeq (,$(filter USE_SFR,$(CONFIGVARS)))
$(error EEOS_SF, AGORA_SF, and JEANS_SF all require USE_SFR)
endif
endif

# Enforce only one SF model at a time
SF_MODELS := $(filter EEOS_SF AGORA_SF JEANS_SF,$(CONFIGVARS))
ifneq ($(word 2,$(SF_MODELS)),)
$(error Only one SF model may be active at a time. Currently enabled: $(SF_MODELS))
endif

ifneq (,$(filter EEOS_SF,$(CONFIGVARS)))
OBJS  += star_formation/sfr_eEOS.o
endif

ifneq (,$(filter AGORA_SF,$(CONFIGVARS)))
OBJS  += star_formation/sfr_AGORA.o
endif

ifneq (,$(filter JEANS_SF,$(CONFIGVARS)))
OBJS  += star_formation/sfr_JEANS.o
endif
endif

#INDIVIDUAL_STAR_BY_STAR_FORMATION
ifneq (,$(filter INDIVIDUAL_STAR_BY_STAR_FORMATION,$(CONFIGVARS)))
ifeq (,$(filter USE_SFR,$(CONFIGVARS)))
    $(error INDIVIDUAL_STAR_BY_STAR_FORMATION requires USE_SFR)
endif
ifeq (,$(filter STAR_PARTICLES 2,$(CONFIGVARS)))
    $(error INDIVIDUAL_STAR_BY_STAR_FORMATION requires STAR_PARTICLES=2)
endif

OBJS += star_formation/individual_star_formation/sfr_starbystar.o \
        star_formation/individual_star_formation/individual_star_formation.o \
        star_formation/individual_star_formation/sf_starbystar.o \
        star_formation/individual_star_formation/sf_massdrain.o
SUBDIRS += star_formation/individual_star_formation
endif

#STARS
ifneq (,$(filter STARS,$(CONFIGVARS)))
OBJS += stars/star.o 
INCL += stars/star.h
SUBDIRS += stars 
endif

define add_define
  grep -qxF '#define $(1)' $(BUILD_DIR)/arepoconfig.h || echo '#define $(1)' >> $(BUILD_DIR)/arepoconfig.h
endef

ifneq (,$(filter STAR_FEEDBACK,$(CONFIGVARS)))
CONFIGVARS += WINDS RADIATION SUPERNOVAE
$(shell $(call add_define,WINDS))
$(shell $(call add_define,SUPERNOVAE))
endif

ifneq (,$(filter RADIATION,$(CONFIGVARS)))
CONFIGVARS += RADIATION_PRESSURE PHOTOELECTRIC_HEATING DISSOCIATION PHOTOIONIZATION
$(shell $(call add_define,RADIATION_PRESSURE))
$(shell $(call add_define,PHOTOELECTRIC_HEATING))
$(shell $(call add_define,DISSOCIATION))
$(shell $(call add_define,PHOTOIONIZATION))
endif

STAR_FEEDBACK_ACTIVE = WINDS RADIATION_PRESSURE PHOTOELECTRIC_HEATING DISSOCIATION PHOTOIONIZATION SUPERNOVAE

STAR_RADIATION_ACTIVE = RADIATION_PRESSURE PHOTOELECTRIC_HEATING DISSOCIATION PHOTOIONIZATION  

ifneq (,$(filter $(STAR_FEEDBACK_ACTIVE),$(CONFIGVARS)))
CONFIGVARS += STAR_FEEDBACK_ACTIVE
$(shell $(call add_define,STAR_FEEDBACK_ACTIVE))
endif

ifneq (,$(filter $(STAR_RADIATION_ACTIVE),$(CONFIGVARS)))
CONFIGVARS += STAR_RADIATION_ACTIVE
$(shell $(call add_define,STAR_RADIATION_ACTIVE))
endif

ifneq (,$(filter STAR_FEEDBACK_ACTIVE,$(CONFIGVARS)))
ifeq (,$(filter STAR_PARTICLES,$(CONFIGVARS)))
$(error STAR_FEEDBACK_ACTIVE requires STAR_PARTICLES)
endif
endif

ifneq (,$(filter STAR_PARTICLES STAR_FEEDBACK_ACTIVE,$(CONFIGVARS)))
ifeq (,$(filter STARS,$(CONFIGVARS)))
$(error STAR_PARTICLES or STAR_FEEDBACK_ACTIVE requires STARS)
endif
endif

ifneq (,$(filter STAR_PARTICLES,$(CONFIGVARS)))
ifeq (,$(filter USE_SFR STAR_FEEDBACK_ACTIVE,$(CONFIGVARS)))
$(error STAR_PARTICLES requires USE_SFR or STAR_FEEDBACK_ACTIVE)
endif
endif

ifneq (,$(filter STAR_HOST_REFINEMENT,$(CONFIGVARS)))
ifeq (,$(filter STAR_FEEDBACK_ACTIVE,$(CONFIGVARS)))
$(error STAR_HOST_REFINEMENT requires STAR_FEEDBACK_ACTIVE)
endif
endif

ifneq (,$(filter STAR_RADIATION_ACTIVE,$(CONFIGVARS)))
ifeq ($(strip $(findstring GRACKLE_CHEMISTRY 2,$(CONFIGVARS))$(findstring GRACKLE_CHEMISTRY 3,$(CONFIGVARS))),)
$(error STAR_RADIATION_ACTIVE requires GRACKLE_CHEMISTRY >= 2)
endif
endif

ifneq (,$(filter STAR_PARTICLES,$(CONFIGVARS)))
OBJS += stars/star_particle.o
INCL += stars/star_particle.h  
endif

ifneq (,$(filter STAR_FEEDBACK_ACTIVE,$(CONFIGVARS)))
OBJS += stars/star_density.o \
        stars/star_update.o \
        stars/star_interpolation.o \
        stars/star_tables.o
INCL += stars/star_proto.h \
        stars/star_tables.h
endif

ifneq (,$(filter WINDS SUPERNOVAE,$(CONFIGVARS)))
OBJS += stars/star_feedback.o 
endif

ifneq (,$(filter STAR_RADIATION_ACTIVE,$(CONFIGVARS)))
OBJS += extern/chealpix.o \
        stars/star_radiation.o \
        stars/star_radiation_tree.o
INCL += extern/chealpix.h \
        stars/star_radiation.h
SUBDIRS += extern
endif

ifneq (,$(filter SUPERNOVAE,$(CONFIGVARS)))
ifeq (,$(filter TREE_BASED_TIMESTEPS,$(CONFIGVARS)))
$(warning SUPERNOVAE without TREE_BASED_TIMESTEPS does not limit timesteps)
endif
endif

#BLACKHOLES
ifneq (,$(filter BLACKHOLES,$(CONFIGVARS)))
OBJS += blackholes/bh.o
INCL += blackholes/bh.h
SUBDIRS += blackholes
endif

ifneq (,$(filter BH_FEEDBACK,$(CONFIGVARS)))
CONFIGVARS += BH_THERMAL_FEEDBACK BH_JET_FEEDBACK 
endif

BH_ACCRETION_ACTIVE = BONDI_ACCRETION TORQUE_ACCRETION ADP_ACCRETION

BH_FEEDBACK_ACTIVE = BH_THERMAL_FEEDBACK BH_JET_FEEDBACK 

ifneq (,$(filter $(BH_ACCRETION_ACTIVE),$(CONFIGVARS)))
CONFIGVARS += BH_ACCRETION_ACTIVE
$(shell $(call add_define,BH_ACCRETION_ACTIVE))
endif

ifneq (,$(filter $(BH_FEEDBACK_ACTIVE),$(CONFIGVARS)))
CONFIGVARS += BH_FEEDBACK_ACTIVE
$(shell $(call add_define,BH_FEEDBACK_ACTIVE))
endif

ifneq (,$(filter BH_ACCRETION_ACTIVE BH_FEEDBACK_ACTIVE,$(CONFIGVARS)))
CONFIGVARS += BH_ACTIVE
$(shell $(call add_define,BH_ACTIVE))
endif

ifneq (,$(filter BH_ACTIVE,$(CONFIGVARS)))
ifeq (,$(filter BLACKHOLES,$(CONFIGVARS)))
$(error BH_ACTIVE requires BLACKHOLES)
endif
endif

ifneq (,$(filter BH_JET_FEEDBACK,$(CONFIGVARS)))
ifneq (,$(filter BH_CONSTANT_RADIUS,$(CONFIGVARS)))
$(error BH_CONSTANT_RADIUS DOES NOT WORK WITH BH_JET_FEEDBACK YET)
endif
endif

# Ensure at most one BH accretion model is enabled
ACCR_SELECTED := $(filter $(BH_ACCRETION_ACTIVE),$(CONFIGVARS))
ACCR_SELECTED_COUNT := $(words $(ACCR_SELECTED))
ifneq (,$(shell test $(ACCR_SELECTED_COUNT) -gt 1 2>/dev/null && echo yes))
$(error ONLY ONE ACCRETION MODEL MAY BE ACTIVE; FOUND: $(ACCR_SELECTED))
endif

ifneq (,$(filter BH_ACTIVE,$(CONFIGVARS)))
OBJS += blackholes/bh_density.o \
        blackholes/bh_update.o 
INCL += blackholes/bh_proto.h 
endif

ifneq (,$(filter BH_ACCRETION_ACTIVE,$(CONFIGVARS)))
OBJS += blackholes/bh_accretion.o \
        blackholes/bh_swallow.o   
endif

ifneq (,$(filter BH_FEEDBACK_ACTIVE,$(CONFIGVARS)))
OBJS += blackholes/bh_feedback.o   
endif

ifeq (FOF,$(findstring FOF,$(CONFIGVARS)))
OBJS += fof/fof.o \
        fof/fof_distribute.o \
        fof/fof_findgroups.o \
        fof/fof_io.o \
        fof/fof_nearest.o \
        fof/fof_sort_kernels.o \
        fof/fof_vars.o
INCL += fof/fof.h
SUBDIRS += fof

ifeq (FIND_HALOS,$(findstring FIND_HALOS,$(CONFIGVARS)))
OBJS += fof/fof_seeding.o 
endif
endif

endif

ifeq (HALO_SEEDING,$(findstring HALO_SEEDING,$(CONFIGVARS)))
OBJS    += fof/fof_seeding.o \
		   fof/fof_seeding_registry.o
INCL    += fof/fof_seeding.h
endif
endif

ifeq (SUBFIND,$(findstring SUBFIND,$(CONFIGVARS)))
OBJS += subfind/subfind.o \
        subfind/subfind_vars.o \
        subfind/subfind_serial.o \
        subfind/subfind_coll_tree.o \
        subfind/subfind_properties.o \
        subfind/subfind_so.o \
        subfind/subfind_distribute.o \
        subfind/subfind_collective.o \
        subfind/subfind_findlinkngb.o \
        subfind/subfind_nearesttwo.o \
        subfind/subfind_loctree.o \
        subfind/subfind_coll_domain.o \
        subfind/subfind_coll_treewalk.o \
        subfind/subfind_density.o \
        subfind/subfind_io.o \
        subfind/subfind_sort_kernels.o \
        subfind/subfind_reprocess.o \
        subfind/subfind_so_potegy.o
INCL += subfind/subfind.h
SUBDIRS += subfind
endif

ifeq (BLACKHOLE_SEEDING,$(findstring BLACKHOLE_SEEDING,$(CONFIGVARS)))
OBJS    += blackholes/bh_seed.o
INCL    += blackholes/bh_proto.h
SUBDIRS += blackholes
endif

##########################
#combine compiler options#
##########################

CFLAGS = $(OPTIMIZE) $(MPICH_INCL) $(HDF5_INCL) -DH5_USE_16_API $(GSL_INCL) $(GMP_INCL) $(FFTW_INCL) $(HWLOC_INCL) $(CELIB_INCL) $(GRACKLE_INCL) -I$(BUILD_DIR)

LIBS = $(GMP_LIB) $(MPICH_LIB) $(HDF5_LIB) -lhdf5 -lz $(GSL_LIB) $(FFTW_LIB) $(HWLOC_LIB) $(CELIB_LIB) $(MATH_LIB) $(GRACKLE_LIB) $(LDFLAGS)

FOPTIONS = $(OPTIMIZE)
FFLAGS = $(FOPTIONS)


SUBDIRS := $(addprefix $(BUILD_DIR)/,$(SUBDIRS))
OBJS := $(addprefix $(BUILD_DIR)/,$(OBJS)) $(BUILD_DIR)/compile_time_info.o $(BUILD_DIR)/compile_time_info_hdf5.o $(BUILD_DIR)/version.o
INCL := $(addprefix $(SRC_DIR)/,$(INCL)) $(BUILD_DIR)/arepoconfig.h

TO_CHECK := $(addsuffix .check, $(OBJS) $(patsubst $(SRC_DIR)%, $(BUILD_DIR)%, $(INCL)) )
TO_CHECK +=  $(BUILD_DIR)/Makefile.check
CONFIG_CHECK = $(BUILD_DIR)/$(notdir $(CONFIG)).check

DOCS_CHECK = $(BUILD_DIR)/README.check

################
#create subdirs#
################
RESULT := $(shell mkdir -p $(SUBDIRS))


#############
#build rules#
#############

all: check build

build: $(EXEC)

$(EXEC): $(OBJS)
	$(LINKER) $(OPTIMIZE) $(OBJS) $(LIBS) -o $(EXEC) 

lib$(LIBRARY).a: $(filter-out $(BUILD_DIR)/main/main.o,$(OBJS))
	$(AR) -rcs lib$(LIBRARY).a $(OBJS)

clean:
	@echo Cleaning all build files...
	@rm -f $(OBJS) $(EXEC) lib$(LIBRARY).a
	@rm -f $(BUILD_DIR)/compile_time_info.c $(BUILD_DIR)/compile_time_info_hdf5.c $(BUILD_DIR)/arepoconfig.h
	@rm -f $(BUILD_DIR)/version.c
	@rm -f $(TO_CHECK) $(CONFIG_CHECK)
	@rm -rf $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(INCL) $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/compile_time_info.o: $(BUILD_DIR)/compile_time_info.c $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/compile_time_info_hdf5.o: $(BUILD_DIR)/compile_time_info_hdf5.c $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu $(INCL) $(MAKEFILES)
	$(NVCC)  -c $< -o $@

# sanity checks:

check: $(CONFIG_CHECK)

check_docs: $(DOCS_CHECK)

$(CONFIG_CHECK): $(TO_CHECK) $(CONFIG) check.py
	@$(PYTHON) check.py 2 $(CONFIG) $(CONFIG_CHECK) defines_extra $(TO_CHECK)

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.c Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.F
	touch $@

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.f90
	touch $@

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.F90
	touch $@

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.cc
	touch $@

$(BUILD_DIR)/%.h.check: $(SRC_DIR)/%.h Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.o.check: $(BUILD_DIR)/%.c Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.h.check: $(BUILD_DIR)/%.h Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/Makefile.check: Makefile Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 3 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/Config.check: Template-Config.sh check.py
	@$(PYTHON) check.py 4 Template-Config.sh $@