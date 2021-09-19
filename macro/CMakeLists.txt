# CMakeLists.txt for FAS package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(FemtoDstAnalyzer)

# You need to tell CMake where to find the ROOT installation. This can be done
# in a number of ways:
#   - ROOT built with classic configure/make use the provided
#   $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix
#   for ROOT
if (ROOT_CMAKE)
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
else (ROOT_CMAKE)
  set(ROOT_PREINSTALLED "/usr/lib64/Geant4-10.3.0/Modules")
  list(APPEND CMAKE_MODULE_PATH ${ROOT_PREINSTALLED})
endif (ROOT_CMAKE)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
if (ROOT_CMAKE)
  find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net
	  HINTS ENV{PATH})
else (ROOT_CMAKE)
  find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
endif (ROOT_CMAKE)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
if (ROOT_CMAKE)
  include(${ROOT_USE_FILE})
endif (ROOT_CMAKE)


#---Define _VANILLA_ROOT_ variable needed for the project
#add_definitions(-D_VANILLA_ROOT_)

set(ST_FEMTO_DST_LIB libStFemtoDst.so)
set(ST_FEMTO_DST_CUSTOM_DIR $ENV{ST_FEMTO_DST_BUILD_DIR})
set(ST_FEMTO_DST_CUSTOM_INC_DIR $ENV{ST_FEMTO_DST_INC_DIR})

find_library(ST_FEMTO_DST REQUIRED
	NAMES ${ST_FEMTO_DST_LIB}
	HINTS ${ST_FEMTO_DST_CUSTOM_DIR})

set(INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}
  ${ST_FEMTO_DST_CUSTOM_INC_DIR}
  ${ROOT_INCLUDE_DIRS}
)

include_directories(${INCLUDE_DIRECTORIES})
set (CMAKE_CXX_STANDARD 11)
add_definitions(${ROOT_CXX_FLAGS})

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")

#---Make executable

add_executable(FemtoDstAnalyzer_Hadrons main_Hadrons.cpp)
target_link_libraries(FemtoDstAnalyzer_Hadrons ${ST_FEMTO_DST} ${ROOT_LIBRARIES})
add_executable(FemtoDstAnalyzer_PID main_PID.cpp)
target_link_libraries(FemtoDstAnalyzer_PID ${ST_FEMTO_DST} ${ROOT_LIBRARIES})
add_executable(PIDcomb main_PIDcomb.cpp)
target_link_libraries(PIDcomb ${ST_FEMTO_DST} ${ROOT_LIBRARIES})
add_executable(combPID_Gaus main_FitCombPID.cpp)
target_link_libraries(combPID_Gaus ${ST_FEMTO_DST} ${ROOT_LIBRARIES})