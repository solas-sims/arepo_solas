# CMakeLists.txt to build the AREPO library/binaries

# SUBDIRS = . \
#           debug_md5 \
#           domain \
#           gitversion \
#           gravity \
#           gravity/pm \
#           hydro \
#           init \
#           io \
#           main \
#           mesh \
#           mesh/voronoi \
#           mpi_utils \
#           ngbtree \
#           pm \
#           star_formation \
#           time_integration \
#           utils \
set(AREPO_SOURCES
    debug_md5/calc_checksum.c
    debug_md5/Md5.c
    domain/domain.c
    domain/domain_balance.c
    domain/domain_box.c
    domain/domain_counttogo.c
    domain/domain_DC_update.c
    domain/domain_exchange.c
    domain/domain_rearrange.c
    domain/domain_sort_kernels.c
    domain/domain_toplevel.c
    domain/domain_vars.c
    domain/peano.c
    # gitversion/version.c
    gravity/accel.c
    gravity/forcetree.c
    gravity/forcetree_ewald.c
    gravity/forcetree_optimizebalance.c
    gravity/forcetree_walk.c
    gravity/grav_external.c
    gravity/grav_softening.c
    gravity/gravdirect.c
    gravity/gravtree.c
    gravity/gravtree_forcetest.c
    gravity/longrange.c
    gravity/pm/pm_periodic2d.c
    gravity/pm/pm_periodic.c
    gravity/pm/pm_mpi_fft.c
    gravity/pm/pm_nonperiodic.c
    hydro/finite_volume_solver.c
    hydro/gradients.c
    hydro/riemann.c
    hydro/riemann_hllc.c
    hydro/riemann_hlld.c
    hydro/scalars.c
    hydro/update_primitive_variables.c
    init/begrun.c
    init/density.c
    init/init.c
    io/global.c
    io/hdf5_util.c
    io/io.c
    io/io_fields.c
    io/logs.c
    io/parameters.c
    io/read_ic.c
    io/restart.c
    main/allvars.c
    main/main.c
    main/run.c
    main/update.c
    mesh/criterion_derefinement.c
    mesh/criterion_refinement.c
    mesh/refinement.c
    mesh/set_vertex_velocities.c
    mesh/voronoi/voronoi.c
    mesh/voronoi/voronoi_1d.c
    mesh/voronoi/voronoi_1d_spherical.c
    mesh/voronoi/voronoi_3d.c
    mesh/voronoi/voronoi_check.c
    mesh/voronoi/voronoi_derefinement.c
    mesh/voronoi/voronoi_dynamic_update.c
    mesh/voronoi/voronoi_exchange.c
    mesh/voronoi/voronoi_ghost_search.c
    mesh/voronoi/voronoi_gradients_lsf.c
    mesh/voronoi/voronoi_gradients_onedims.c
    mesh/voronoi/voronoi_refinement.c
    mesh/voronoi/voronoi_utils.c
    mpi_utils/checksummed_sendrecv.c
    mpi_utils/hypercube_allgatherv.c
    mpi_utils/mpi_util.c
    mpi_utils/myalltoall.c
    mpi_utils/sizelimited_sendrecv.c
    mpi_utils/pinning.c
    ngbtree/ngbtree.c
    ngbtree/ngbtree_search.c
    ngbtree/ngbtree_walk.c
    star_formation/sfr_eEOS.c
    star_formation/starformation.c
    star_formation/sfr_AGORA.c
    star_formation/sfr_JEANS.c
    time_integration/darkenergy.c
    time_integration/do_gravity_hydro.c
    time_integration/driftfac.c
    time_integration/predict.c
    time_integration/timestep.c
    time_integration/timestep_treebased.c
    utils/allocate.c
    utils/debug.c
    utils/mpz_extension.c
    utils/mymalloc.c
    utils/parallel_sort.c
    utils/system.c
)

set(AREPO_INCLUDES_DIRS
    debug_md5
    domain
    gitversion
    gravity
    gravity/pm
    main
    mesh
    mesh/voronoi
    time_integration
    utils
    mpi_utils
    ngbtree
)

if (AREPO_ENABLE_TWODIMS)
    list(APPEND AREPO_SOURCES
        mesh/voronoi/voronoi_2d.c
    )
endif()

if (AREPO_ENABLE_MPI)
    list(APPEND AREPO_SOURCES
        mpi_utils/myIBarrier.c
    )
    list(APPEND AREPO_INCLUDES_DIRS
        mpi_utils
    )
endif()

if (AREPO_ENABLE_MHD)
    list(APPEND AREPO_SOURCES
        hydro/mhd.c
    )
    list(APPEND AREPO_INCLUDES_DIRS
        hydro
    )
endif()

if (AREPO_ENABLE_ADDBACKGROUNDGRID)
    list(APPEND AREPO_SOURCES
        add_backgroundgrid/add_bggrid.c
        add_backgroundgrid/calc_weights.c
        add_backgroundgrid/distribute.c
    )
    list(APPEND AREPO_INCLUDES_DIRS
        add_backgroundgrid
    )

endif()

if (AREPO_ENABLE_COOLING)
    list(APPEND AREPO_SOURCES
        cooling/cooling.c
    )
    list(APPEND AREPO_INCLUDES_DIRS
        cooling
    )
endif()

if (AREPO_ENABLE_USE_GRACKLE)
    list(APPEND AREPO_SOURCES
        cooling/grackle.c
    )
    list(APPEND AREPO_INCLUDES_DIRS
        cooling
    )
endif()

if (AREPO_ENABLE_FOF)
    list(APPEND AREPO_SOURCES
        fof/fof.c
        fof/fof_distribute.c
        fof/fof_findgroups.c
        fof/fof_io.c
        fof/fof_nearest.c
        fof/fof_sort_kernels.c
        fof/fof_vars.c
    )
    list(APPEND AREPO_INCLUDES_DIRS
        fof
    )
endif()

if (AREPO_ENABLE_SUBFIND)
    list(APPEND AREPO_SOURCES
        subfind/subfind.c
        subfind/subfind_vars.c
        subfind/subfind_serial.c
        subfind/subfind_coll_tree.c
        subfind/subfind_properties.c
        subfind/subfind_so.c
        subfind/subfind_distribute.c
        subfind/subfind_collective.c
        subfind/subfind_findlinkngb.c
        subfind/subfind_nearesttwo.c
        subfind/subfind_loctree.c
        subfind/subfind_coll_domain.c
        subfind/subfind_coll_treewalk.c
        subfind/subfind_density.c
        subfind/subfind_io.c
        subfind/subfind_sort_kernels.c
        subfind/subfind_reprocess.c
        subfind/subfind_so_potegy.c
    )
    list(APPEND AREPO_INCLUDES_DIRS
        subfind
    )
endif()

if (AREPO_ENABLE_BLACKHOLES)
    list(APPEND AREPO_SOURCES
        blackholes/bh_density.c
        blackholes/bh_jet_density.c
        blackholes/bh_feedback.c
        blackholes/bh_refinement.c
    )
    list(APPEND AREPO_INCLUDES_DIRS
        blackholes
    )
endif() 

if (AREPO_ENABLE_STARS)
    list(APPEND AREPO_SOURCES
        stars/star_density.c
        stars/star_feedback.c 
    )
    list(APPEND AREPO_INCLUDES_DIRS
        stars
    )
endif()

if (AREPO_ENABLE_USE_CELIB)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/celib/src ${CMAKE_CURRENT_SOURCE_DIR}/celib/data)
endif()

if (AREPO_ENABLE_PROFILE_UTIL)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../profile_util/include)
endif()

list(APPEND AREPO_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/../compile_time_info.c
    ${CMAKE_CURRENT_BINARY_DIR}/../compile_time_info_hdf5.c
    ${CMAKE_CURRENT_BINARY_DIR}/../gitversion.c
)

# ifeq (MYIBARRIER,$(findstring MYIBARRIER,$(CONFIGVARS)))
# OBJS    += mpi_utils/myIBarrier.o
# INCL    += mpi_utils/myIBarrier.h
# endif

# ifeq (MHD,$(findstring MHD,$(CONFIGVARS)))
# OBJS    += hydro/mhd.o
# endif

# ifeq (ADDBACKGROUNDGRID,$(findstring ADDBACKGROUNDGRID,$(CONFIGVARS)))
# OBJS    += add_backgroundgrid/add_bggrid.o \
#            add_backgroundgrid/calc_weights.o \
#            add_backgroundgrid/distribute.o
# INCL    += add_backgroundgrid/add_bggrid.h
# SUBDIRS += add_backgroundgrid
# endif

# ifeq (COOLING,$(findstring COOLING,$(CONFIGVARS)))
# OBJS    += cooling/cooling.o
# INCL    += cooling/cooling_vars.h \
#            cooling/cooling_proto.h
# SUBDIRS += cooling
# endif

# ifeq (USE_GRACKLE,$(findstring USE_GRACKLE,$(CONFIGVARS)))
# OBJS    += cooling/grackle.o
# endif

# ifeq (FOF,$(findstring FOF,$(CONFIGVARS)))
# OBJS    += fof/fof.o \
#            fof/fof_distribute.o \
#            fof/fof_findgroups.o \
#            fof/fof_io.o \
#            fof/fof_nearest.o \
#            fof/fof_sort_kernels.o \
#            fof/fof_vars.o
# INCL    += fof/fof.h
# SUBDIRS += fof
# endif

# ifeq (SUBFIND,$(findstring SUBFIND,$(CONFIGVARS)))
# OBJS	+= subfind/subfind.o \
#            subfind/subfind_vars.o \
#            subfind/subfind_serial.o \
#            subfind/subfind_coll_tree.o \
#            subfind/subfind_properties.o \
#            subfind/subfind_so.o \
#            subfind/subfind_distribute.o \
#            subfind/subfind_collective.o \
#            subfind/subfind_findlinkngb.o \
#            subfind/subfind_nearesttwo.o \
#            subfind/subfind_loctree.o \
#            subfind/subfind_coll_domain.o \
#            subfind/subfind_coll_treewalk.o \
#            subfind/subfind_density.o \
#            subfind/subfind_io.o \
#            subfind/subfind_sort_kernels.o \
#            subfind/subfind_reprocess.o \
#            subfind/subfind_so_potegy.o
# INCL	+= subfind/subfind.h
# SUBDIRS += subfind
# endif

# ifeq (BLACKHOLES,$(findstring BLACKHOLES,$(CONFIGVARS)))
# OBJS    += blackholes/bh_density.o \
#            blackholes/bh_jet_density.o\
#            blackholes/bh_feedback.o \
#            blackholes/bh_refinement.o
# INCL    += blackholes/bh_proto.h
# SUBDIRS += blackholes
# endif

# ifeq (STARS,$(findstring STARS,$(CONFIGVARS)))
# OBJS    += stars/star_density.o \
#            stars/star_feedback.o 
# INCL    += stars/star_proto.h
# SUBDIRS += stars
# endif

# ???

# INCL += debug_md5/Md5.h \
#         domain/bsd_tree.h \
#         domain/domain.h \
#         gitversion/version.h\
#         gravity/forcetree.h \
#         main/allvars.h \
#         main/proto.h \
#         mesh/mesh.h \
#         mesh/voronoi/voronoi.h \
#         time_integration/timestep.h \
#         utils/dtypes.h \
#         utils/generic_comm_helpers2.h \
#         utils/timer.h
# )

# # add default subdirectories
# add_subdirectory(debug_md5)
# add_subdirectory(domain)
# add_subdirectory(gitversion)
# add_subdirectory(gravity)
# add_subdirectory(gravity)
# add_subdirectory(hydro)
# add_subdirectory(init)
# add_subdirectory(io)
# add_subdirectory(main)
# add_subdirectory(mesh)
# add_subdirectory(ngbtree)
# add_subdirectory(pm)
# add_subdirectory(star_formation)
# add_subdirectory(time_integration)
# add_subdirectory(utils)


include_directories(${AREPO_INCLUDES_DIRS})
add_library(arepolib STATIC ${AREPO_SOURCES})
target_compile_definitions(arepolib PUBLIC ${AREPO_DEFINES})

# if (AREPO_CXX_FLAGS)
# 	target_compile_options(arepolib PUBLIC ${AREPO_CXX_FLAGS})
# endif()
set_target_properties(arepolib PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
sanitizer_options(arepolib)


add_executable(arepo main/main.c)
target_link_libraries(arepo arepolib ${AREPO_LIBS})
if (AREPO_ENABLE_USE_CELIB)
    target_link_libraries(arepo celib)
endif()
if (AREPO_ENABLE_PROFILE_UTIL)
    target_link_libraries(arepo profile_util)
    add_compile_definitions(arepo PUBLIC "-DENABLE_PROFILE_UTIL")
endif()
if (AREPO_LINK_FLAGS)
	set_target_properties(arepo PROPERTIES LINK_FLAGS ${AREPO_LINK_FLAGS})
endif()
if (AREPO_ENABLE_FFTW)
    if (AREPO_DOUBLEPRECISION_FFTW)
       target_link_libraries(arepo fftw3 fftw3_mpi)
    else()
         target_link_libraries(arepo fftw3f fftw3f_mpi)
    endif()
endif()
set_target_properties(arepo PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
target_compile_definitions(arepo PUBLIC ${AREPO_DEFINES})

if (BUILD_TESTING)
	add_subdirectory(tests)
endif()
