cmake_minimum_required(VERSION 2.6)

project(mediocr-dicom)

# Our current version
# When doing backwards compatible changes, raise the MINOR version
# When doing backwards incompatible changes, raise MAJOR and set MINOR to 0
set(mediocr-dicom_VERSION_MAJOR 1)
set(mediocr-dicom_VERSION_MINOR 5)

# Use our own FindDCMTK.cmake, in which library linking order is fixed
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})

# Include DCMTK
find_package(DCMTK REQUIRED)

add_library(mediocr-dicom src/dicom_file.cpp src/dicom_coders.cpp)
target_link_libraries(mediocr-dicom ${DCMTK_LIBRARIES})
target_compile_options(mediocr-dicom PUBLIC "-std=c++11")
target_compile_definitions(mediocr-dicom PUBLIC "-DHAVE_CONFIG_H") # Fix DCMTK Compilation
target_include_directories(mediocr-dicom SYSTEM PUBLIC ${DCMTK_INCLUDE_DIRS})
target_include_directories(mediocr-dicom
	INTERFACE
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<INSTALL_INTERFACE:"..">
)

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
