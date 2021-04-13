# Include external modules
include(GNUInstallDirs)
include(ExternalProject)

# If `psi` is not found by `find_package`
if(NOT TARGET psi::psi)
  message(STATUS "Using bundled psi library")
  set(psi_SOURCE_DIR ${PROJECT_SOURCE_DIR}/psi)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${psi_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	#set(BUILD_WITH_BDSG OFF CACHE BOOL "Build with BDSG support")
	set(BUILD_PSIKT OFF CACHE BOOL "Build PSI seeding/indexing command-line tool (`psi kt`)")
	set(USE_BUNDLED_ALL ON CACHE BOOL "Use all bundled dependencies")
  add_subdirectory(${psi_SOURCE_DIR} EXCLUDE_FROM_ALL)
endif()

if(ParallelHashmap_FOUND)
  set(ParallelHashmap_PKG_CFLAGS "-I${ParallelHashmap_INCLUDE_DIRS}")
endif(ParallelHashmap_FOUND)

# If `parallel-hashmap` is not found by `find_package`
if(NOT TARGET ParallelHashmap::ParallelHashmap)
  message(STATUS "Using bundled parallel-hashmap library")
  set(ParallelHashmap_SOURCE_DIR ${PROJECT_SOURCE_DIR}/parallel-hashmap)
  ExternalProject_Add(parallelhashmap_git
    DOWNLOAD_COMMAND git -C ${PROJECT_SOURCE_DIR} submodule update --init --recursive -- ${ParallelHashmap_SOURCE_DIR}
    SOURCE_DIR ${ParallelHashmap_SOURCE_DIR}
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
  ExternalProject_Get_Property(parallelhashmap_git INSTALL_DIR)
  add_library(ParallelHashmap::ParallelHashmap INTERFACE IMPORTED)
  set_target_properties(ParallelHashmap::ParallelHashmap PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${INSTALL_DIR}/include>;$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>")
endif()

# If `concurrentqueue` is not found by `find_package`
if(NOT TARGET concurrentqueue)
  message(STATUS "Using bundled concurrentqueue library")
  set(CONCURRENTQUEUE_SOURCE_DIR ${PROJECT_SOURCE_DIR}/concurrentqueue)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${CONCURRENTQUEUE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(${CONCURRENTQUEUE_SOURCE_DIR} EXCLUDE_FROM_ALL)
endif()

# If `BBHash` is not found by `find_package`
if(NOT TARGET BBHash::BBHash)
  message(STATUS "Using bundled BBHash library")
  set(BBHash_SOURCE_DIR ${PROJECT_SOURCE_DIR}/BBHash)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${BBHash_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_library(BBHash::BBHash INTERFACE IMPORTED)
  set_target_properties(BBHash::BBHash PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${BBHash_SOURCE_DIR}")
endif()

# If `zstr` is not found by `find_package`
if(NOT TARGET zstr::zstr)
  message(STATUS "Using bundled zstr library")
  set(ZSTR_SOURCE_DIR ${PROJECT_SOURCE_DIR}/zstr)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${ZSTR_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_library(zstr::zstr INTERFACE IMPORTED)
  set_target_properties(zstr::zstr PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${ZSTR_SOURCE_DIR}/src")
endif()
