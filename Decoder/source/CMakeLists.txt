# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_VERBOSE_MAKEFILE 1)

project(ROOT_COMPILE)

set(Targets NI2Root Debug NIEventShow NIAnal NIVis)
#set(Targets ToRootFour Show EventShow Offset Height PeakTime RiseFallTime)
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)


#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED)
#find_package(ROOT REQUIRED COMPONENTS RIO Net)
include_directories(${ROOT_INCLUDE_DIRS})


#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)


file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*)
file(GLOB headers ${PROJECT_SOURCE_DIR}/Event/*)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*)


#---Create a shared library with geneated dictionary
foreach(path IN LISTS Targets)
  get_filename_component(target ${path} NAME_WE)
  add_executable(${target} ${path}
													 ${sources} ${headers})
	target_link_libraries(${target} 
												${ROOT_LIBRARIES}
												Spectrum
												libstdc++.so
												libm.so.6)
endforeach()


set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/..)
install(TARGETS ${Targets} DESTINATION bin)
