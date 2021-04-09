# Find 'jemalloc' library.
#
# This set the following variables:
#   - jemalloc_FOUND
#   - jemalloc_VERSION
#   - jemalloc_INCLUDE_DIRS
#   - jemalloc_LIBRARIES
#
# and the following imported targets:
#   - jemalloc::jemalloc

if(jemalloc_INCLUDE_DIRS)
  set(jemalloc_FIND_QUIETLY TRUE)
else()
  # Try pkg-config, first.
  find_package(PkgConfig QUIET)
  pkg_check_modules(jemalloc QUIET jemalloc>=5.0.1)
endif(jemalloc_INCLUDE_DIRS)

# If jemalloc_INCLUDE_DIRS is still not set, run `jemalloc-config` to get the config.
if(NOT jemalloc_INCLUDE_DIRS)
  execute_process(COMMAND jemalloc-config --includedir OUTPUT_VARIABLE jemalloc_INCLUDE_DIRS)
  execute_process(COMMAND jemalloc-config --libs OUTPUT_VARIABLE jemalloc_LIBRARIES)
  execute_process(COMMAND jemalloc-config --libdir OUTPUT_VARIABLE jemalloc_LIBRARY_DIRS)
  execute_process(COMMAND jemalloc-config --version OUTPUT_VARIABLE jemalloc_VERSION)
endif(NOT jemalloc_INCLUDE_DIRS)

## handle the QUIETLY and REQUIRED arguments and set jemalloc_FOUND to TRUE if
## all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(jemalloc DEFAULT_MSG jemalloc_INCLUDE_DIRS jemalloc_LIBRARIES jemalloc_LIBRARY_DIRS)

mark_as_advanced(jemalloc_FOUND jemalloc_VERSION jemalloc_INCLUDE_DIRS jemalloc_LIBRARIES jemalloc_LIBRARY_DIRS)

# Define `jemalloc::jemalloc` imported target
if(jemalloc_FOUND AND NOT TARGET jemalloc::jemalloc)
  add_library(jemalloc::jemalloc INTERFACE IMPORTED)
  set_target_properties(jemalloc::jemalloc PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${jemalloc_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${jemalloc_LIBRARIES}"
    INTERFACE_LINK_DIRECTORIES "${jemalloc_LIBRARY_DIRS}")
endif()
