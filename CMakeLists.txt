cmake_minimum_required(VERSION 3.15..3.27)

project(adaptive_remesh VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

set(CMAKE_BUILD_TYPE "Release")  # change to Debug for debugging
#set(CMAKE_BUILD_TYPE "Debug")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(EXT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext)

#add_compile_definitions(IGL_PARALLEL_FOR_FORCE_SERIAL)


if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()


set(Python_FIND_VIRTUALENV FIRST)
find_package(Python 3.8 COMPONENTS Interpreter ${DEV_MODULE} NumPy REQUIRED)


# Detect the installed nanobind package and import it into CMake
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)


# putting this before fetching igl makes igl's cmake scripts able to see this
# Eigen instead of trying to fetch its own Eigen
find_package(Threads REQUIRED)
find_package(Eigen3 3.4...<4.0)
# this find_package(Eigen3) is not REQUIRED but if it passes (found Eigen >=3.4), 
# igl won't try to fetch its own Eigen. But if this doesn't pass, (either no eigen or
# older version), fetching igl will fetch its own eigen and make it available for us


# fetch libigl
include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.5.0
    SOURCE_DIR ${EXT_DIR}/libigl
)
FetchContent_MakeAvailable(libigl)

include_directories(
	# libigl
	${libigl_SOURCE_DIR}/include
	)
	
add_library(adaptive_remesh_lib
	SHARED

	# Headers
	src/collapse_edges.h
	src/equalize_valences.h
	src/adaptive_remesh_botsch.h
	src/split_edges.h
	src/split_edges_until_bound.h
	src/tangential_smoothing.h

	# Source
	# src/collapse_edges.cpp
	# src/equalize_valences.cpp
	src/adaptive_remesh_botsch.cpp
	# src/split_edges.cpp
	# src/split_edges_until_bound.cpp
	# src/tangential_smoothing.cpp
	)

target_link_libraries(adaptive_remesh_lib Threads::Threads Eigen3::Eigen)
target_compile_options(adaptive_remesh_lib PRIVATE -O3)


nanobind_add_module(bkremeshlerps src/nanobind_bindings.cpp)
target_include_directories(bkremeshlerps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(bkremeshlerps PUBLIC adaptive_remesh_lib)


