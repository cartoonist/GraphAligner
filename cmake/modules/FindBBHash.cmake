# Find 'BBHash' library.
#
# This set the following variables:
#   - BBHash_FOUND
#   - BBHash_INCLUDE_DIRS
#
# and the following imported targets:
#   - BBHash::BBHash

if(BBHash_INCLUDE_DIRS)
  set(BBHash_FIND_QUIETLY TRUE)
else()
  # If BBHash_INCLUDE_DIRS is not set, this searches for the header file.
  find_path(BBHash_INCLUDE_DIRS BooPHF.h)
endif(BBHash_INCLUDE_DIRS)

# handle the QUIETLY and REQUIRED arguments and set BBHash_FOUND to
# TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BBHash DEFAULT_MSG BBHash_INCLUDE_DIRS)

mark_as_advanced(BBHash_INCLUDE_DIRS)

# Define `BBHash::BBHash` imported target
if(BBHash_FOUND AND NOT TARGET BBHash::BBHash)
  add_library(BBHash::BBHash INTERFACE IMPORTED)
  set_target_properties(BBHash::BBHash PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${BBHash_INCLUDE_DIRS}")
endif()
