cmake_minimum_required(VERSION 2.6)

project(mediocr-dicom)

# make {#include <mediocr-dicom/...>} work inside the repository
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# C++11 mode
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# Use our own FindDCMTK.cmake, in which library linking order is fixed
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})

# Include DCMTK
find_package(DCMTK REQUIRED)
include_directories(SYSTEM ${DCMTK_INCLUDE_DIRS})
# Fix DCMTK compilation
add_definitions("-DHAVE_CONFIG_H")

add_library(mediocr-dicom src/dicom_file.cpp src/dicom_coders.cpp)
target_link_libraries(mediocr-dicom ${DCMTK_LIBRARIES})

install(TARGETS mediocr-dicom EXPORT mediocr-dicom-targets ARCHIVE DESTINATION lib)
install(FILES
	mediocr-dicom/dicom_file.hpp
	mediocr-dicom/dicom_coders.hpp
	DESTINATION "include/mediocr-dicom")
