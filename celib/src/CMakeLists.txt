# CMakeLists.txt to build the celib
set(CELIB_SOURCES
    IMFTypes.c 
    Info.c 
    WriteIMF.c 
    Init.c 
    YieldsInfo.c 
    LifeTime.c 
    IMFFunctions.c 
    Integral.c 
    SNIIYields.c 
    SNIaYields.c 
    SNIIRate.c 
    SNIaRate.c 
    NSMYields.c 
    NSMRate.c 
    InitialMetallicity.c 
    AGBMassLoss.c 
    UnifiedAPIs.c 
    RunParameters.c 
    CheckRunParameters.c 
    Utilities.c
)

set(CELIB_INCLUDES_DIRS
    src
)

# set_property(SOURCE ${CELIB_SOURCES} PROPERTY LANGUAGE CXX)
add_library(celib STATIC ${CELIB_SOURCES} )
# target_compile_features(celib PRIVATE cxx_std_11)
target_include_directories(celib PUBLIC ${CELIB_INCLUDES_DIRS})
set_target_properties(celib PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
