# Find 'Mummer' library.
#
# This set the following variables:
#   - mummer_FOUND
#   - mummer_VERSION
#   - mummer_INCLUDE_DIRS
#   - mummer_LIBRARIES
#
# and the following imported targets:
#   - mummer::mummer

if(mummer_INCLUDE_DIRS)
  set(mummer_FIND_QUIETLY TRUE)
else()
  # Try pkg-config, first.
  find_package(PkgConfig QUIET)
  pkg_check_modules(mummer QUIET mummer>=4.0.0)
  # If mummer_INCLUDE_DIRS is not set, this searches for the header/library file.
  find_path(mummer_INCLUDE_DIRS mummer/sparseSA.hpp)
  find_library(mummer_LIBRARIES umdmummer)
endif(mummer_INCLUDE_DIRS)

## handle the QUIETLY and REQUIRED arguments and set mummer_FOUND to TRUE if
## all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(mummer DEFAULT_MSG mummer_INCLUDE_DIRS mummer_LIBRARIES)

mark_as_advanced(mummer_FOUND mummer_VERSION mummer_INCLUDE_DIRS mummer_LIBRARIES)

# Define `mummer::mummer` imported target
if(mummer_FOUND AND NOT TARGET mummer::mummer)
  add_library(mummer::mummer INTERFACE IMPORTED)
  set_target_properties(mummer::mummer PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${mummer_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${mummer_LIBRARIES}")
endif()
