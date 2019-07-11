cmake_minimum_required(VERSION 2.6)
project (emConst)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

#Non ho la definizione nel .cpp, non lo metto!
#add_library(sum_alg sum_alg)

add_library(gamma_const gamma_const.cpp)
target_link_libraries(gamma_const  ${ROOT_LIBRARIES})

add_executable(un_emconst un_emconst.cpp)
target_link_libraries(un_emconst gamma_const ${ROOT_LIBRARIES})

