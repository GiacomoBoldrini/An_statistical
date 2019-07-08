cmake_minimum_required(VERSION 2.6)
project (pt_pl)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

add_library(monte_carlo monte_carlo.cpp)
target_link_libraries(monte_carlo ${ROOT_LIBRARIES})

add_executable(main main.cpp)
target_link_libraries(main monte_carlo ${ROOT_LIBRARIES})