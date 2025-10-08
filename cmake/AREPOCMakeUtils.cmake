#
# How we find MPI and set it up
#
macro(arepo_find_mpi)
    message("MPI enabled, finding MPI ... ")
    find_package(MPI REQUIRED)
    if (MPI_FOUND)
        include_directories(${MPI_CXX_INCLUDE_PATH})
        include_directories(${MPI_C_INCLUDE_PATH})
        link_libraries(${MPI_CXX_LIBRARIES})
        link_libraries(${MPI_C_LIBRARIES})
        list(APPEND CMAKE_CXX_FLAGS ${MPI_CXX_FLAGS})
        list(APPEND CMAKE_C_FLAGS ${MPI_C_FLAGS})
        list(APPEND CMAKE_DEFINES "_MPI")
        set(AREPO_HAS_MPI "Yes")
    else()
        message(SEND_ERROR "MPI enabled but not found. Please check configuration or disable MPI")
    endif()
endmacro()

macro(arepo_mpi)
    set(AREPO_HAS_MPI No)
    if (AREPO_ENABLE_MPI)
        arepo_find_mpi()
    endif()
endmacro()

macro(arepo_find_openmp)
    message("OpenMP enabled, finding OpenMP ... ")
    find_package(OpenMP)
    #FindOpenMP()
    if (OpenMP_FOUND)
        # this is a silly hack for CUDA, just sets the openmp flag to -fopenmp
        if (ENABLE_OPENMP_COMPILE_FLAG_HACK)
            set(OpenMP_CXX_FLAGS "-fopenmp")
        endif()
        set(AREPO_HAS_OPENMP "Yes")
        add_link_options(${OpenMP_CXX_FLAGS})
    else()
        message(SEND_ERROR "OpenMP enabled but not found. Please check configuration or disable OpenMP")
    endif()
endmacro()

macro(arepo_openmp)
    set(AREPO_HAS_OPENMP No)
    if (AREPO_ENABLE_OPENMP)
        arepo_find_openmp()
    endif()
endmacro()

#
# How we find HDF5 and set it up
#
macro(arepo_find_hdf5)
	# FindHDF5 needs an environment variable, oddly, unlike
	# most other packages that use normal cmake variables
	if (HDF5_ROOT)
		set(ENV{HDF5_ROOT} ${HDF5_ROOT})
	endif()
	find_package(HDF5 COMPONENTS C REQUIRED)
	if (HDF5_FOUND)
    	include_directories(${HDF5_INCLUDE_DIRS})
	    link_libraries(${HDF5_LIBRARIES})
	    # target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES})
		set(AREPO_HAS_HDF5 Yes)
		#check if parallel hdf present
		if (HDF5_IS_PARALLEL AND AREPO_HAS_MPI AND AREPO_ALLOWPARALLELHDF5)
			set (ENV{HDF5_PREFER_PARALLEL} true)
			set(AREPO_HAS_PARALLEL_HDF5 Yes)
			list(APPEND AREPO_DEFINES USEPARALLELHDF)
			if (HDF5_VERSION VERSION_GREATER "1.10.0" AND AREPO_ALLOWCOMPRESSIONPARALLELHDF5)
				set(AREPO_HAS_COMPRESSED_HDF5 Yes)
				list(APPEND AREPO_DEFINES USEHDFCOMPRESSION)
				list(APPEND AREPO_DEFINES PARALLELCOMPRESSIONACTIVE)
			endif()
		else()
			if (AREPO_ALLOWCOMPRESSIONHDF5)
				set(AREPO_HAS_COMPRESSED_HDF5 Yes)
				list(APPEND AREPO_DEFINES USEHDFCOMPRESSION)
			endif()
		endif()
    endif()
endmacro()

macro(arepo_hdf5)
    set(AREPO_HAS_HDF5 No)
    set(AREPO_HAS_COMPRESSED_HDF5 No)
    set(AREPO_HAS_PARALLEL_HDF5 No)
    if (AREPO_ENABLE_HDF5)
    	arepo_find_hdf5()
    endif()
endmacro()

macro(arepo_gsl)
	find_package(GSL REQUIRED)
	include_directories(${GSL_INCLUDE_DIRS})
    link_libraries(${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})
	#target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})
endmacro()

macro(arepo_gmp)
    if (DEFINED ENV{GMP_ROOT})
        set(GMP_ROOT $ENV{GMP_ROOT})
    else()
        set(GMP_ROOT "/usr/local" CACHE PATH "Path to GMP installation")
    endif()
    message(STATUS "Using GMP_ROOT: ${GMP_ROOT}")
    if (NOT EXISTS "${GMP_ROOT}/include/gmp.h")
        message(FATAL_ERROR "GMP not found in GMP_ROOT: ${GMP_ROOT}. Please set the GMP_ROOT environment variable or CMake cache entry to point to your GMP installation.")
    endif()
    include_directories(${GMP_ROOT}/include)
    link_directories(${GMP_ROOT}/lib)
	# target_include_directories(${GMP_ROOT}/include)
	link_libraries(gmp)
	# target_link_libraries(${PROJECT_NAME} gmp)
endmacro()

# this might require an update to deal with cufft and hipfft
macro(arepo_find_fftw)
	find_package(FFTW REQUIRED)
	include_directories(${FFTW_INCLUDE_DIRS})
	link_libraries(${FFTW_LIBRARIES})
	# target_link_libraries(${PROJECT_NAME} ${FFTW_LIBRARIES})
endmacro()

macro(arepo_fftw)
    if (AREPO_ENABLE_FFTW)
    	arepo_find_fftw()
    endif()
endmacro()

macro(arepo_find_hip)
    # Find hip
    message("HIP enabled, finding HIP ... ")
    find_package(HIP)
    if (HIP_FOUND)
        # Link with HIP
        enable_language(HIP)
        if(NOT DEFINED CMAKE_HIP_STANDARD)
            set(CMAKE_HIP_STANDARD 17)
            set(CMAKE_HIP_STANDARD_REQUIRED ON)
        endif()
        list(APPEND AREPO_DEFINES "_HIP")
        set(AREPO_HAS_HIP "Yes")
        add_compile_options("-fPIE")
        add_compile_options("-D_HIP")
        target_link_libraries(hip::device)
        set(CMAKE_HIP_FLAGS ${CMAKE_CXX_FLAGS} ${CMAKE_HIP_FLAGS})
        if (ENABLE_HIP_AMD)
            message(STATUS "Using AMD HIP ")
            add_definitions("-D__HIP_PLATFORM_AMD__")
            set(GPU_TARGETS "gfx90a" CACHE STRING "GPU targets to compile for")
        endif()
    else()
        message(SEND_ERROR "HIP enabled but not found. Please check configuration or disable HIP")
    endif()
endmacro()

macro(arepo_find_cuda)
    # Find cuda
    message(STATUS "CUDA enabled, finding CUDA ... ")
    find_package(CUDA REQUIRED)
    if (CUDA_FOUND)
        enable_language(CUDA)
        list(APPEND AREPO_DEFINES "_CUDA")
        add_definitions("-D_CUDA")
        set(AREPO_HAS_CUDA "Yes")
        message(STATUS "SYSTEM_ARCH: ${SYSTEM_ARCH}")
        if(SYSTEM_ARCH STREQUAL "aarch64")
            set(CUDA_ARCH_DIR "sbsa-linux")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -O3 -ftree-vectorize")
        else()
            set(CUDA_ARCH_DIR "x86_64-linux")
        endif()

        if(NOT DEFINED CMAKE_CUDA_STANDARD)
            set(CMAKE_CUDA_STANDARD 17)
            set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        endif()

        # ------------------------------------------------------------------------------
        # CUDA compile options and Thrust settings
        # ------------------------------------------------------------------------------
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --expt-extended-lambda")
        add_definitions(-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA)
        add_definitions(-DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_CPP)

    else()
        message(SEND_ERROR "CUDA enabled but not found. Please check configuration or disable CUDA")
    endif()
endmacro()


macro(arepo_hip)
    set(AREPO_HAS_HIP No)
    if (AREPO_ENABLE_HIP)
        arepo_find_hip()
    endif()
endmacro()

macro(arepo_cuda)
    set(AREPO_HAS_CUDA No)
    if (AREPO_ENABLE_CUDA)
        arepo_find_cuda()
    endif()
endmacro()


macro(arepo_report feature)

    # Output feature name and underscore it in the next line
    message("\n${feature}")
    string(REGEX REPLACE "." "-" _underscores ${feature})
    message("${_underscores}\n")

    set(_args "${ARGN}")
    # message(STATUS "_args: ${_args}")
    list(LENGTH _args _nargs)
    math(EXPR _nargs "${_nargs} - 1")
    foreach(_idx RANGE 0 ${_nargs} 2)
        # Items in the list come with a message first, then the variable name
        list(GET _args ${_idx} _msg)
        math(EXPR _idx2 "${_idx} + 1")
        list(GET _args ${_idx2} _varname)
        # message(STATUS "_varname: ${_varname}")

        # We try to keep things up to 80 cols
        string(LENGTH ${_msg} _len)
        math(EXPR _nspaces "75 - ${_len}")
        string(RANDOM LENGTH ${_nspaces} _spaces)
        string(REGEX REPLACE "." " " _spaces "${_spaces}")
        string(CONCAT _msg "${_msg}" ${_spaces})
        message(" ${_msg} ${AREPO_${_varname}}")
    endforeach()
endmacro()

macro(arepo_compilation_summary)
    message(STATUS "AREPO successfully configured with the following settings:")
    arepo_report("MPI-specifics"
            "MPI support" ENABLE_MPI
            )
    arepo_report("File formats"
            "HDF5" ENABLE_HDF5
            "Compressed HDF5" HAS_COMPRESSED_HDF5
            "Parallel HDF5" HAS_PARALLEL_HDF5
    )
    arepo_report("OpenMP-specifics"
            "OpenMP support" ENABLE_OPENMP
            "OpenMP GPU offloading support" ENABLE_OPENMPGPU
            )
    arepo_report("GPU-specifics"
            "CUDA support" ENABLE_CUDA
            "HIP support" ENABLE_HIP)
    arepo_report("Particle-specifics"
            "Use stars" ENABLE_STARS
            "Use black holes" ENABLE_BLACKHOLES
            )
    arepo_report("Hydro-specifics"
            "No Hydrodynamics" ENABLE_NOHYDRO
            "Magnetohydrodynamics (MHD)" ENABLE_MHD
            "Riemann solver: HLLC" ENABLE_RIEMANN_HLLC
            "Riemann solver: HLLD" ENABLE_RIEMANN_HLLD
        )

	message("")
	message("Compilation")
	message("=========================")
	message("Compiler ID          : ${CMAKE_CXX_COMPILER_ID}" )
	message("Compiler             : ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ID}" )
	message("Build Type           : ${CMAKE_BUILD_TYPE}")
	message("Build Flags          : ${${CMAKE_BUILD_FLAGS}}")
	if(NOT ENABLE_SANITIZER STREQUAL "none")
		message("Sanitizer            : -fsanitize=${ENABLE_SANITIZER}")
	endif()
	message("----------------------")
	message("Include directories  : ${AREPO_INCLUDE_DIRS}")
	message("AREPO Macros         : ${AREPO_DEFINES}")
	message("AREPO Lib            : ${AREPO_LIBS}")
	message("AREPO C++ flags      : ${AREPO_C_FLAGS}")
	message("AREPO Link flags     : ${AREPO_LINK_FLAGS}")
	message("")
endmacro()

#
# Some useful utils
#

#
# list values as bullet points
#
function(list_to_bulletpoints result)
    list(REMOVE_AT ARGV 0)
    set(temp "")
    foreach(item ${ARGV})
        set(temp "${temp}* ${item}")
    endforeach()
    set(${result} "${temp}" PARENT_SCOPE)
endfunction(list_to_bulletpoints)

#
# valid the option choosen based on allowed values
#
function(validate_option name values)
    string(TOLOWER ${${name}} needle_lower)
    string(TOUPPER ${${name}} needle_upper)
    list(FIND ${values} ${needle_lower} IDX_LOWER)
    list(FIND ${values} ${needle_upper} IDX_UPPER)
    if(${IDX_LOWER} LESS 0 AND ${IDX_UPPER} LESS 0)
        list_to_bulletpoints(POSSIBLE_VALUE_LIST ${${values}})
        message(FATAL_ERROR "\n########################################################################\n"
                            "Invalid value '${${name}}' for option ${name}\n"
                            "Possible values are : "
                            "${POSSIBLE_VALUE_LIST}"
                            "\n"
                            "########################################################################")
    endif()
endfunction(validate_option)

#
# Function to add a "doc" target, which will doxygen out the given
# list of directories
#
function(try_add_doc_target doc_dirs)
	find_program(DOXYGEN_FOUND doxygen)
	if (NOT DOXYGEN_FOUND)
		return()
	endif()

	# Create a target for each individual doc directory, then a final one
	# that depends on them all
	message("-- Adding doc target for directories: ${doc_dirs}")
	set(_dependencies "")
	set(x 1)
	foreach(_doc_dir IN ITEMS ${doc_dirs})
		add_custom_command(OUTPUT ${_doc_dir}/xml/index.xml
		                   COMMAND doxygen
		                   WORKING_DIRECTORY ${_doc_dir})
		add_custom_target(doc_${x}
		                  COMMAND doxygen
		                  WORKING_DIRECTORY ${_doc_dir})
		list(APPEND _dependencies doc_${x})
		math(EXPR x "${x} + 1")
	endforeach()
	add_custom_target(doc DEPENDS "${_dependencies}")
endfunction()

#
# - Prevent in-source builds.
# https://stackoverflow.com/questions/1208681/with-cmake-how-would-you-disable-in-source-builds/
#
macro(prevent_in_source_builds)
    # make sure the user doesn't play dirty with symlinks
    get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
    get_filename_component(srcdir2 "${CMAKE_SOURCE_DIR}/.." REALPATH)
    get_filename_component(srcdir3 "${CMAKE_SOURCE_DIR}/../src" REALPATH)
    get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

    # disallow in-source builds
    if("${srcdir}" STREQUAL "${bindir}" OR "${srcdir2}" STREQUAL "${bindir}" OR "${srcdir3}" STREQUAL "${bindir}")
        message(FATAL_ERROR "\
            CMake must not to be run in the source directory. 
            Rather create a dedicated build directory and run CMake there. 
            To clean up after this aborted in-place compilation:
            rm -r CMakeCache.txt CMakeFiles
        ")
    endif()
endmacro()

#
# set the default build and also store the compilation flags
# as a string based on the currently choosen flags
#
macro(my_set_build_type)
	set(default_build_type "Release")
	if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
		message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
		set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
			STRING "Choose the type of build." FORCE)
		# Set the possible values of build type for cmake-gui
		set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
			"Debug" "Release" "MinSizeRel" "RelWithDebInfo")
	endif()
	#set(ACTIVE_COMPILE_OPTIONS )
endmacro()

macro(enable_santizer_option)
    set(ENABLE_SANITIZER "none" CACHE STRING "Select a code sanitizer option (none (default), address, leak, thread, undefined)")
    mark_as_advanced(ENABLE_SANITIZER)
    set(ENABLE_SANITIZER_VALUES none address leak thread undefined)
    set_property(CACHE ENABLE_SANITIZER PROPERTY STRINGS ${ENABLE_SANITIZER_VALUES})
    validate_option(ENABLE_SANITIZER ENABLE_SANITIZER_VALUES)
    string(TOLOWER ${ENABLE_SANITIZER} ENABLE_SANITIZER)
endmacro()

function(sanitizer_options mytarget)
    if(NOT ENABLE_SANITIZER STREQUAL "none")
        if((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
            target_compile_options(${mytarget} PUBLIC -fsanitize=${ENABLE_SANITIZER})
            target_link_options(${mytarget} PUBLIC -fsanitize=${ENABLE_SANITIZER})
        else()
            message(WARNING "ENABLE_SANITIZER option not supported by ${CMAKE_CXX_COMPILER_ID} compilers. Ignoring.")
            set(ENABLE_SANITIZER "none")
        endif()
    endif()
endfunction()

# ... (Include necessary modules)
include(FindPythonInterp REQUIRED)

# Custom function to read YAML file and set variables
function(read_yaml FILE)
  find_package(PythonInterp REQUIRED)
  message(STATUS "Reading ${FILE}")
  execute_process(
    COMMAND ${CMAKE_CURRENT_LIST_DIR}/cmake/parse_yaml.py ${FILE} ${CMAKE_CURRENT_BINARY_DIR} 
  )
endfunction()


macro(arepo_configure_from_yaml)
    if (NOT EXISTS "${AREPO_YAML_CONFIG}")
        message(FATAL_ERROR "Configuration file ${AREPO_YAML_CONFIG} not found.")
    endif()
    # Read config.yaml and set variables
    read_yaml("${AREPO_YAML_CONFIG}")
    include(${CMAKE_CURRENT_BINARY_DIR}/generated_options.cmake)
endmacro()


#run some macros automatically
prevent_in_source_builds()
my_set_build_type()
enable_santizer_option()
