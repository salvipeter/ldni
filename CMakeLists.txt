cmake_minimum_required(VERSION 3.1.0)

project(ldni)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
endif()

set(LIBGEOM_DIR "../libgeom") # https://github.com/salvipeter/libgeom
include_directories("${LIBGEOM_DIR}" ${EIGEN3_INCLUDE_DIR})
link_directories("${LIBGEOM_DIR}/release")

add_library(ldni STATIC ldni.cc cell-point.cc)

target_link_libraries(ldni "geom")

add_executable(ldni-test ldni-test.cc)

target_link_libraries(ldni-test ldni)
