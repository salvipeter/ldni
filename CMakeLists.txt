cmake_minimum_required(VERSION 3.1.0)

project(ldni)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif(CMAKE_COMPILER_IS_GNUCXX)

set(Transfinite_DIR "../transfinite/src") # https://github.com/salvipeter/transfinite
include_directories("${Transfinite_DIR}/geom" ${EIGEN3_INCLUDE_DIR})
link_directories("${Transfinite_DIR}/geom/Release" "${Transfinite_DIR}/../release/geom")

add_library(ldni STATIC ldni.cc cell-point.cc)

target_link_libraries(ldni "geom")

add_executable(ldni-test ldni-test.cc)

target_link_libraries(ldni-test ldni)
