cmake_minimum_required(VERSION 2.6)

project(mediocr-dicom)

# Our current version
# When doing backwards compatible changes, raise the MINOR version
# When doing backwards incompatible changes, raise MAJOR and set MINOR to 0
set(mediocr-dicom_VERSION_MAJOR 1)
set(mediocr-dicom_VERSION_MINOR 2)

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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/mediocr-dicom-config.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/mediocr-dicom-config.cmake @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/mediocr-dicom-config-version.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/mediocr-dicom-config-version.cmake @ONLY)

install(TARGETS mediocr-dicom EXPORT mediocr-dicom-targets ARCHIVE DESTINATION lib)
install(FILES
	mediocr-dicom/dicom_file.hpp
	mediocr-dicom/dicom_coders.hpp
	DESTINATION "include/mediocr-dicom")
install(FILES
	cmake/Usemediocr-dicom.cmake
	${CMAKE_CURRENT_BINARY_DIR}/mediocr-dicom-config.cmake
	${CMAKE_CURRENT_BINARY_DIR}/mediocr-dicom-config-version.cmake
	cmake/modules/FindPackageHandleStandardArgs.cmake
	DESTINATION lib/cmake/mediocr-dicom/)
install(EXPORT mediocr-dicom-targets DESTINATION lib/cmake/mediocr-dicom/)
