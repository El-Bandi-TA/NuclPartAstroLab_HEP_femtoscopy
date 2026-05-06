cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

project(HEP_Femtoscopy)

add_definitions(-DHEPQC_ROOT="${CMAKE_SOURCE_DIR}")
add_definitions(-DHEP_DATA="/project/femtoscopy/data/PHENIX_AuAu_200GeV")

find_package(ROOT REQUIRED)
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

add_library(QuantumCorrLib SHARED ${LIB_SOURCES})
target_link_libraries(QuantumCorrLib ${ROOT_LIBRARIES})

add_executable(analyzeCzphi apps/analyzeCzphi.cpp)
target_link_libraries(analyzeCzphi QuantumCorrLib ${ROOT_LIBRARIES})

add_executable(analyzeCQ apps/analyzeCQ.cpp)
target_link_libraries(analyzeCQ QuantumCorrLib ${ROOT_LIBRARIES})
