cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

project(HEP_Femtoscopy)

add_definitions(-DHEPQC_ROOT="${CMAKE_SOURCE_DIR}")
add_definitions(-DHEP_DATA="/project/femtoscopy/data/PHENIX_AuAu_200GeV")

find_package(ROOT REQUIRED COMPONENTS Minuit2 MathCore)
find_package(Boost REQUIRED)
execute_process(COMMAND root-config --cxx --version OUTPUT_VARIABLE ROOT_CXX_VERSION)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR}/include)

# file(GLOB LIB_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)
set(
    LIB_SOURCES
    src/event.cpp
    src/kt_pair_cuts.cpp
    src/kt_parser.cpp
    src/particle_tree.cpp
    src/pion.cpp
    src/pool.cpp
    src/utils.cpp
)

file(
    GLOB COULCORR_SOURCES
    ${CMAKE_SOURCE_DIR}/external/CoulCorrLevyIntegral/src/*.cpp
)

add_library(QuantumCorrLib SHARED ${LIB_SOURCES})
target_link_libraries(QuantumCorrLib ${ROOT_LIBRARIES})

add_library(CoulCorrLevyIntegral_static STATIC ${COULCORR_SOURCES})
target_include_directories(CoulCorrLevyIntegral_static PUBLIC
    "external/CoulCorrLevyIntegral/include"
    ${Boost_INCLUDE_DIRS}
)
target_link_libraries(CoulCorrLevyIntegral_static PUBLIC
    ${ROOT_LIBRARIES}
    ROOT::Minuit2
    Boost::boost
    # ${Boost_LIBRARIES}
)
add_library(CoulCorrLevyIntegral SHARED ${COULCORR_SOURCES})
target_include_directories(CoulCorrLevyIntegral PUBLIC
    "external/CoulCorrLevyIntegral/include"
    ${Boost_INCLUDE_DIRS}
)
target_link_libraries(CoulCorrLevyIntegral PUBLIC
    ${ROOT_LIBRARIES}
    ROOT::Minuit2
    Boost::boost
    # ${Boost_LIBRARIES}
)

set_target_properties(CoulCorrLevyIntegral_static
    PROPERTIES OUTPUT_NAME CoulCorrLevyIntegral
)

add_executable(analyzeCzphi apps/analyzeCzphi.cpp)
target_link_libraries(analyzeCzphi QuantumCorrLib ${ROOT_LIBRARIES})

add_executable(analyzeCQ apps/analyzeCQ.cpp)
target_link_libraries(analyzeCQ QuantumCorrLib ${ROOT_LIBRARIES})

add_executable(fitCQ apps/fitCQ.cpp)
target_link_libraries(fitCQ
    PRIVATE
    QuantumCorrLib
    CoulCorrLevyIntegral
    ${ROOT_LIBRARIES}
    ROOT::Minuit2
    ROOT::MathCore
)