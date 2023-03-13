cmake_minimum_required(VERSION 2.8)

project(./Calibration)

aux_source_directory(./src SRCS)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
message(STATUS "OpenCV version: ${OpenCV_VERSION}")
message(STATUS "include path: ${OpenCV_INCLUDE_DIRS}")

add_executable(test ${SRCS})

target_link_libraries(test ${OpenCV_LIBS} -lpthread)

