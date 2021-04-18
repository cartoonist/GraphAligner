cmake_minimum_required(VERSION 3.10)
project(GraphAligner VERSION 1.1.0 LANGUAGES CXX)

# Options
option(PSI_STATS "Build with PSI statistics" OFF)

# Include external modules
include(GNUInstallDirs)
# Adding CMake module path.
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
# Finding dependencies.
find_package(ZLIB REQUIRED)
find_package(jemalloc REQUIRED)
find_package(sdsl REQUIRED)
find_package(Mummer REQUIRED)
find_package(VGio REQUIRED)
find_package(Boost 1.67 REQUIRED COMPONENTS serialization program_options)
find_package(OpenMP)
# Bundled dependencies
find_package(BBHash QUIET)
#find_package(zstr QUIET)             # No find module, always use bundled
#find_package(ParallelHashmap QUIET)  # Handled by PSI
# FIXME: PSI configure file complains that dependencies targets cannot be found
#find_package(psi QUIET)              # Always use bundled
# Include external dependencies
include(${PROJECT_SOURCE_DIR}/external.cmake)

# Load `CheckGitRevision` module
include(CheckGitRevision)
set(VERSION "${GIT_REVISION} ${GIT_COMMIT_DATE}")

# Defining GraphAligner source files
set(GRAPHALIGNER_SOURCES "src/AlignerMain.cpp"
	"src/Aligner.cpp"
	"src/fastqloader.cpp"
	"src/BigraphToDigraph.cpp"
	"src/ThreadReadAssertion.cpp"
	"src/AlignmentGraph.cpp"
	"src/CommonUtils.cpp"
	"src/GraphAlignerWrapper.cpp"
	"src/GfaGraph.cpp"
	"src/AlignmentCorrectnessEstimation.cpp"
	"src/MummerSeeder.cpp"
	"src/ReadCorrection.cpp"
	"src/MinimizerSeeder.cpp"
	"src/AlignmentSelection.cpp")

# Defining target 'GraphAligner'
add_executable(GraphAligner ${GRAPHALIGNER_SOURCES})
if(TARGET parallelhashmap_git)
	add_dependencies(GraphAligner parallelhashmap_git)
endif()
target_compile_definitions(GraphAligner PRIVATE VERSION="${VERSION}")
target_compile_options(GraphAligner PRIVATE -g -Wall -Wpedantic -Wextra)
# Use C++17
target_compile_features(GraphAligner PRIVATE cxx_std_17)
# Defining include directories
target_include_directories(GraphAligner
  PRIVATE jemalloc::jemalloc
	PRIVATE Mummer::Mummer
	PRIVATE VGio::VGio
	PRIVATE Boost::serialization
	PRIVATE Boost::program_options
# PRIVATE ParallelHashmap::ParallelHashmap
	PRIVATE concurrentqueue
	PRIVATE BBHash::BBHash
	PRIVATE zstr::zstr
  PRIVATE psi::psi)
# Defining link libraries
target_link_libraries(GraphAligner
  PRIVATE jemalloc::jemalloc
	PRIVATE Mummer::Mummer
	PRIVATE VGio::VGio
	PRIVATE Boost::serialization
	PRIVATE Boost::program_options
# PRIVATE ParallelHashmap::ParallelHashmap
	PRIVATE BBHash::BBHash
	PRIVATE zstr::zstr
	PRIVATE ZLIB::ZLIB
  PRIVATE psi::psi)
# Install GraphAligner
install(TARGETS GraphAligner DESTINATION ${CMAKE_INSTALL_BINDIR})

# Defining UntipRelative source files
set(UNTIPRELATIVE_SOURCES "src/UntipRelative.cpp"
	"src/CommonUtils.cpp"
	"src/GfaGraph.cpp"
	"src/fastqloader.cpp"
	"src/ThreadReadAssertion.cpp")

# Defining target 'UntipRelative'
add_executable(UntipRelative ${UNTIPRELATIVE_SOURCES})
target_compile_options(UntipRelative PRIVATE -g -Wall -Wpedantic -Wextra)
# Use C++17
target_compile_features(UntipRelative PRIVATE cxx_std_17)
# Defining include directories
target_include_directories(UntipRelative
	PRIVATE VGio::VGio
	PRIVATE zstr::zstr)
# Defining link libraries
target_link_libraries(UntipRelative
	PRIVATE VGio::VGio
	PRIVATE zstr::zstr
	PRIVATE ZLIB::ZLIB)
# Install UntipRelative
install(TARGETS UntipRelative DESTINATION ${CMAKE_INSTALL_BINDIR})
