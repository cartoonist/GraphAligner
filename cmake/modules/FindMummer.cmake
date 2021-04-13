# Find 'Mummer' library.
#
# This set the following variables:
#   - Mummer_FOUND
#   - Mummer_VERSION
#   - Mummer_INCLUDE_DIRS
#   - Mummer_LIBRARIES
#
# and the following imported targets:
#   - Mummer::Mummer

if(Mummer_INCLUDE_DIRS)
  set(Mummer_FIND_QUIETLY TRUE)
else()
  # Try pkg-config, first.
  find_package(PkgConfig QUIET)
  pkg_check_modules(Mummer QUIET mummer>=4.0.0)
  # If Mummer_INCLUDE_DIRS is not set, this searches for the header/library file.
  find_path(Mummer_INCLUDE_DIRS mummer/sparseSA.hpp)
  find_library(Mummer_LIBRARIES umdmummer)
endif(Mummer_INCLUDE_DIRS)

## handle the QUIETLY and REQUIRED arguments and set Mummer_FOUND to TRUE if
## all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Mummer DEFAULT_MSG Mummer_INCLUDE_DIRS Mummer_LIBRARIES)

mark_as_advanced(Mummer_FOUND Mummer_VERSION Mummer_INCLUDE_DIRS Mummer_LIBRARIES)

# Define `Mummer::Mummer` imported target
if(Mummer_FOUND AND NOT TARGET Mummer::Mummer)
  add_library(Mummer::Mummer INTERFACE IMPORTED)
  set_target_properties(Mummer::Mummer PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${Mummer_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${Mummer_LIBRARIES}")
endif()
