cmake_minimum_required(VERSION 2.6)
project (prng)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

include_directories(/Users/boldrinicoder/RooUnfold/include)
include_directories(/Users/boldrinicoder/RooUnfold/src)
link_directories(/Users/boldrinicoder/RooUnfold/)



add_library(prng prng)

#FOR GSL LIBRARY!
find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(prng GSL::gsl GSL::gslcblas)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)


add_library(xorshiro xorshiro.cpp)
target_link_libraries(xorshiro prng ${ROOT_LIBRARIES})

add_library(xorshiroGen xorshiroGen.cpp)
target_link_libraries(xorshiroGen prng xorshiro ${ROOT_LIBRARIES})

add_library(simulation simulation)
target_link_libraries(simulation prng xorshiro xorshiroGen ${ROOT_LIBRARIES})

add_executable(main main.cpp)
target_link_libraries(main prng xorshiro xorshiroGen simulation ${ROOT_LIBRARIES} RooUnfold Unfold)


#add_executable(un_prngBW_hm un_prngBW_hm.cpp)
# target_link_libraries(un_prngBW_hm prng prngGen)
# 
# add_executable(un_prngBW_inv un_prngBW_inv.cpp)
# target_link_libraries(un_prngBW_inv prng prngGen)
# 
#  add_executable(un_xorshiro128p un_xorshiro128p.cpp)
#  target_link_libraries(un_xorshiro128p prng prngGen xorshiro xorshiroGen)
# 
# 
# add_executable(un_xorshiro128pBWInv un_xorshiro128pBWInv.cpp)
# target_link_libraries(un_xorshiro128pBWInv prng prngGen xorshiro xorshiroGen)

#foreach(EXE un_lgc un_prngBW_hm un_prngBW_inv un_xorshiro128p un_xorshiro128pBWInv)
#    add_executable(${EXE} ${EXE}.cpp)
#    target_link_libraries(prng prngGen xorshiro xorshiroGen)
#endforeach(EXE)

#add_executable(un_lgc un_lgc.cpp)
#target_link_libraries(un_lgc prng prngGen)

#add_executable(un_xorshiro128p un_xorshiro128p.cpp)
#target_link_libraries(un_xorshiro128p prng xorshiro)<
 
