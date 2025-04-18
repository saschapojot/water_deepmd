cmake_minimum_required(VERSION 3.28)
project(water_deepmd)
#set(CMAKE_C_COMPILER "/opt/intel/oneapi/compiler/2025.1/bin/icx")


#set(CMAKE_CXX_COMPILER "/opt/intel/oneapi/compiler/2025.1/bin/icpx")
#set(AVX2_FLAGS "-xCORE-AVX2 -qopt-zmm-usage=high")

message(STATUS "Using C++ compiler: ${CMAKE_CXX_COMPILER}")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CXX_STANDARD_LIB "/usr/lib/x86_64-linux-gnu/libstdc++.so.6")
add_executable(run_mc run_mc.cpp mc_subroutine/mc_read_load_compute.cpp mc_subroutine/mc_read_load_compute.hpp)
target_compile_options(run_mc PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
        $<$<CONFIG:Release>:-ftree-vectorize>
)
set_source_files_properties(
        mc_subroutine/mc_read_load_compute.cpp
        PROPERTIES

        COMPILE_OPTIONS "-O3"
)

set_source_files_properties(
        run_mc.cpp
        PROPERTIES

        COMPILE_OPTIONS "-O3"
)

set(Python3_INCLUDE_DIRS "/home/adada/anaconda3/envs/deepmd_venv/include/python3.10")
set(Python3_LIBRARIES "/home/adada/anaconda3/envs/deepmd_venv/lib/libpython3.10.so")
include_directories(${Python3_INCLUDE_DIRS})
#include_directories("/home/adada/cppLocal/armadillo/include")
include_directories("/home/adada/cppLocal/libtorch_cpu/include")
include_directories("/home/adada/cppLocal/deep_potential_cpp_gcc/include")

#include_directories("/home/adada/cppLocal/deep_potential_cpp_intel/include")

set(LIBTORCH_DIR "/home/adada/cppLocal/libtorch_cpu/lib/libtorch.so")
set(dpcppDir "/home/adada/cppLocal/deep_potential_cpp_gcc/lib/libdeepmd_cc.so")

#set(dpcppDir "/home/adada/cppLocal/deep_potential_cpp_intel/lib/libdeepmd_cc.so")
set(Boost_NO_SYSTEM_PATHS YES)
if (Boost_NO_SYSTEM_PATHS)
    set(Boost_ROOT "/home/adada/cppLocal/boostBuild")
    set(Boost_INCLUDE_DIR "/home/adada/cppLocal/boostBuild/include")
    set(Boost_LIBRARY_DIR "/home/adada/cppLocal/boostBuild/lib")
endif (Boost_NO_SYSTEM_PATHS)

find_package(Boost  REQUIRED COMPONENTS python numpy filesystem)
if(Boost_FOUND)
    message(STATUS "boost root is: ${Boost_ROOT}")
    message(STATUS "boost include path is : ${Boost_INCLUDE_DIR}")
    message(STATUS "boost library path is : ${Boost_LIBRARY_DIR}")
    message(STATUS "boost libraries is : ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIR})
else()
    message(WARNING "boost not found.")
endif()


target_link_libraries(run_mc Boost::filesystem)

target_link_libraries(run_mc Boost::python)
target_link_libraries(run_mc Boost::numpy)
#python
target_link_libraries(run_mc ${Python3_LIBRARIES})

#armadillo
#file(GLOB LIBS_SO "/home/adada/cppLocal/armadillo/lib/libarmadillo.so")
#target_link_libraries(run_mc ${LIBS_SO})

target_link_directories(run_mc PRIVATE  "/home/adada/cppLocal/deep_potential_cpp_gcc/lib")

#target_link_directories(run_mc PRIVATE  "/home/adada/cppLocal/deep_potential_cpp_intel/lib")
target_link_libraries(run_mc ${LIBTORCH_DIR})
target_link_libraries(run_mc ${dpcppDir})

