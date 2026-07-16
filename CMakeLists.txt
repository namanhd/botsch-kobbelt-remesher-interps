cmake_minimum_required(VERSION 3.15..3.27)

project(bkremeshlerps VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)


# Warn if the user invokes CMake directly
if (NOT SKBUILD)
  message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build-core'.
  Running it directly will almost certainly not produce the desired
  result. If you are a user trying to install this package, use the
  command below, which will install all necessary build dependencies,
  compile the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to rerun the above
  after editing C++ files.")
endif()


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
find_package(Python 3.9
  REQUIRED COMPONENTS Interpreter ${DEV_MODULE}
  OPTIONAL_COMPONENTS Development.SABIModule)


# Detect the installed nanobind package and import it into CMake
#execute_process(
#  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
#  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)

find_package(nanobind CONFIG REQUIRED)

find_package(Threads REQUIRED)

# putting this before fetching igl makes igl's cmake scripts able to see this
# Eigen instead of trying to fetch its own Eigen

if(CMAKE_VERSION VERSION_LESS 3.19)
  find_package(Eigen3 3.4...)
else()
  find_package(Eigen3 3.4...<4.0)  # some arbitrarily high version, since 3.4... doesn't work,,
endif()

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


# check if eigen was actually found..
if(Eigen3_FOUND)
    message(STATUS "Eigen3 found successfully!")
    message(STATUS "Eigen3 version: ${Eigen3_VERSION}")
else()
    message(WARNING "Eigen3 was NOT found!")
endif()


include_directories(
	# libigl
	${libigl_SOURCE_DIR}/include
)
	
add_library(bkremeshlerps_lib
	SHARED

	# Headers
	src/collapse_edges.h
	src/equalize_valences.h
	src/adaptive_remesh_botsch.h
	src/split_edges.h
	src/split_edges_until_bound.h
	src/tangential_smoothing.h

	# Source
  # (NOTE the commented out source files are now #included into their headers)
	# src/collapse_edges.cpp
	# src/equalize_valences.cpp
	src/adaptive_remesh_botsch.cpp
	# src/split_edges.cpp
	# src/split_edges_until_bound.cpp
	# src/tangential_smoothing.cpp
)

target_link_libraries(bkremeshlerps_lib Threads::Threads Eigen3::Eigen)
target_compile_options(bkremeshlerps_lib PRIVATE -O3)


nanobind_add_module(
  bkremeshlerps_nb

  # Target the stable ABI for Python 3.12+, which reduces
  # the number of binary wheels that must be built. This
  # does nothing on older Python versions
  STABLE_ABI

  src/nanobind_bindings.cpp
)

target_include_directories(bkremeshlerps_nb PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(bkremeshlerps_nb PUBLIC bkremeshlerps_lib)


# since our code is packaged as a shared library, we have to make the extension module (bkremeshlerps_nb)
# have an RPATH that is "current directory" so that it can find the shared library
if(APPLE)
    set(LINKER_ORIGIN "@loader_path")
else()
    set(LINKER_ORIGIN "$ORIGIN")
endif()

set_target_properties(bkremeshlerps_nb PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${LINKER_ORIGIN}"
)

# Install directive for scikit-build-core
# (have to include the shared library as well)
install(TARGETS bkremeshlerps_lib bkremeshlerps_nb LIBRARY DESTINATION ${PROJECT_NAME})

# Generate type stubs
nanobind_add_stub(
  bkremeshlerps_nb_stub
  MODULE bkremeshlerps_nb
  OUTPUT bkremeshlerps_nb.pyi
  MARKER_FILE py.typed
  PYTHON_PATH $<TARGET_FILE_DIR:bkremeshlerps_nb>
  DEPENDS bkremeshlerps_nb
)

install(FILES ${CMAKE_BINARY_DIR}/bkremeshlerps_nb.pyi ${CMAKE_BINARY_DIR}/py.typed DESTINATION ${PROJECT_NAME})
