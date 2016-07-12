cmake_minimum_required(VERSION 2.8.12)

project(mediocr-dicom)

# Our current version
# When doing backwards compatible changes, raise the MINOR version
# When doing backwards incompatible changes, raise MAJOR and set MINOR to 0
set(mediocr-dicom_VERSION_MAJOR 2)
set(mediocr-dicom_VERSION_MINOR 5)

find_package(wubwubcmake 1.2 REQUIRED)
include(wubwubcmake/sane_install)
include(wubwubcmake/position_independent_code)
include(wubwubcmake/wub_find_modules)
include(wubwubcmake/enable_cpp11)

find_package(DCMTK REQUIRED)

file(GLOB sources "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/mediocr-dicom/*.hpp")
add_library(mediocr-dicom ${sources} ${headers})
target_link_libraries(mediocr-dicom ${DCMTK_LIBRARIES})
target_compile_definitions(mediocr-dicom PUBLIC "-DHAVE_CONFIG_H") # Fix DCMTK Compilation
target_include_directories(mediocr-dicom SYSTEM PUBLIC ${DCMTK_INCLUDE_DIRS})
target_include_directories(mediocr-dicom
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<INSTALL_INTERFACE:.>
)

sane_install(TARGETS mediocr-dicom HEADERS ${headers})
