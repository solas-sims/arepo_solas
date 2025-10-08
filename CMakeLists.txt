cmake_minimum_required(VERSION 3.25)

project(arepo,
  VERSION 0.5
  DESCRIPTION "Utility to report parallel API related information and also timing."
  LANGUAGES C CXX
)

# define cxx 
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)
# add cmake macros
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include (AREPOCMakeUtils)

macro(arepo_option optname optdesc status)
   option(AREPO_${optname} "${optdesc}" "${status}")
endmacro()
arepo_option(ENABLE_MPI "Enable MPI" ON)
arepo_option(ENABLE_OPENMP "Enable OpenMP" ON)
arepo_option(ENABLE_CUDA "Enable CUDA" OFF)
arepo_option(ENABLE_HIP "Enable HIP" OFF)
arepo_option(ENABLE_HIP_AMD "Enable HIP with AMD (ROCM)" ON)
arepo_option(ENABLE_SHARED_LIB "Enable shared library" ON)
arepo_option(ENABLE_TESTS "Enable testing" ON)
arepo_option(ENABLE_CMAKE_PTHREADS_FIX "Enable CMake pthreads fix for ARM" OFF)

set(AREPO_YAML_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/template_config.yaml)

# arepo_option(ENABLE_TWODIMS "Enable 2D simulations" OFF)
# arepo_option(ENABLE_MYIBARRIER "Enable custom MPI barrier" OFF)
# arepo_option(ENABLE_MHD "Enable MHD" OFF)
# arepo_option(ENABLE_ADDBACKGROUNDGRID "Enable adding background grid" OFF)
# arepo_option(ENABLE_COOLING "Enable cooling" OFF)
# arepo_option(ENABLE_USE_GRACKLE "Enable Grackle cooling library" OFF)
# arepo_option(ENABLE_FOF "Enable FOF group finder" OFF)
# arepo_option(ENABLE_SUBFIND "Enable SUBFIND substructure finder" OFF)
# arepo_option(ENABLE_BLACKHOLES "Enable black hole physics" OFF)
# arepo_option(ENABLE_STARS "Enable star formation and feedback" OFF)

include(CTest)
enable_testing()


if ((AREPO_ENABLE_HIP) AND (AREPO_ENABLE_CUDA))
   message(FATAL_ERROR "Build conflict: Requesting both CUDA and HIP. Must have one or the other")
endif()

# check for MPI 
arepo_mpi()
# check for openmp 
arepo_openmp()
# check for hdf5
arepo_hdf5()
# check for gsl 
arepo_gsl()
# check for gmp
arepo_gmp()
# check for fftw
arepo_fftw()
# check for cuda 
arepo_cuda()
# check for hip
arepo_hip()

add_compile_options(-Wall -Wextra -Wpedantic -fPIC)

arepo_configure_from_yaml()

# weird error for cmake on arm it seems 
# see https://github.com/alicevision/geogram/issues/2 
# for similar issue and this resolution
# there might be other reasons to include it. Further testing is required
if (ENABLE_CMAKE_PTHREADS_FIX)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
endif()


arepo_compilation_summary()

#
# Add git revision information in a generated .cpp file, which we then include
# as part of the library
#
# include(cmake/GetGitRevisionDescription.cmake)
# get_git_head_revision(GIT_REFSPEC GIT_SHA1)
# git_local_changes(GIT_HAS_LOCAL_CHANGES)
# set(git_revision_cpp "${CMAKE_CURRENT_BINARY_DIR}/git_revision.cpp")
# configure_file("src/git_revision.cpp.in" "${git_revision_cpp}" @ONLY)

# set the paths for building
# set(AREPO_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(AREPO_SOURCE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
# include_directories(${AREPO_INCLUDE_DIRS})
include_directories(${AREPO_SOURCE_DIRS})
add_subdirectory(${AREPO_SOURCE_DIRS})

# Are we part of an outside build?
# If so, we want to let the caller know what to include
set(_export OFF)
get_directory_property(_hasParent PARENT_DIRECTORY)
if (_hasParent)
	set(_export ON)
endif()

