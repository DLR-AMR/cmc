cmake_minimum_required(VERSION 3.23)

project(cmc
        DESCRIPTION "Lossy and lossless data compression based on AMR"
        LANGUAGES C CXX
        VERSION 0.1.0)

include(GNUInstallDirs)
set(CMC_VERSION "${PROJECT_VERSION}")

# Configuration options
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
option(ENABLE_MPI "Enable features that rely on MPI" ON)
option(WITH_NETCDF "Enable features that rely on netCDF" OFF)
option(WITH_T8CODE "Link against an available t8code installation" ON)
option(ENABLE_ALL_WARNINGS "Use compilation flags -Wall -Wextra -Wpedantic" [OFF])

# Define some paths for potential library linkage
set(CMC_NETCDF_ROOT_DIR "" CACHE PATH "Directory containing netCDF installation")

include(CTest)
enable_testing()

# Check if a shared or a static library should be built
if(BUILD_SHARED_LIBS)
    add_library(cmc SHARED)
else()
    add_library(cmc STATIC)
endif()
add_library(cmc::cmc ALIAS cmc)

# Add some library properties regarding versioning
set_target_properties(cmc PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(cmc PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
set_target_properties(cmc PROPERTIES OUTPUT_NAME cmc)

# Set C and C++ standard for compilations
target_compile_features(cmc PUBLIC c_std_11 cxx_std_20)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set a default build type if none was specified
set(default_build_type "Release")

# Check if a build type has been given, if not we build a release configuration 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build. Build types available: Release Debug RelWithDebInfo" FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()

# Add flags for a release configuration
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set (cmc_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
    set (cmc_CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
endif()

# Add flags for debug configuration
if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    set (cmc_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
    set (cmc_CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
    target_compile_definitions(cmc PUBLIC ENABLE_DEBUG=1)
    set(CMC_ENABLE_DEBUG ON)
endif()

# Check if MPI is available if activated
if(ENABLE_MPI)
    find_package(MPI COMPONENTS C REQUIRED)
    target_compile_definitions(cmc PUBLIC ENABLE_MPI=1)
    target_link_libraries(cmc MPI::MPI_C)
    if(MPI_C_FOUND)
        message("-- An MPI installation has been found")
    endif ()
    set(CMC_ENABLE_MPI ON)
else()
    set(CMC_ENABLE_MPI OFF)
    target_compile_definitions(cmc PUBLIC ENABLE_MPI=0)
endif()

# Check if netCDF is available if activated
if(WITH_NETCDF)
    include(CheckCSourceCompiles)
    # Check if a netCDF library path has been set explicitly
    if (CMC_NETCDF_ROOT_DIR)
        find_library(NETCDF netcdf HINTS ${CMC_NETCDF_ROOT_DIR}/lib REQUIRED)
        message("-- Found netCDF installation: ${NETCDF}")
        set(CMAKE_REQUIRED_INCLUDES "${CMC_NETCDF_ROOT_DIR}/include")
        set(CMAKE_REQUIRED_LIBRARIES "${NETCDF}")
        set(CMAKE_REQUIRED_LINK_DIRECTORIES "${CMC_NETCDF_ROOT_DIR}/lib")
        check_c_source_compiles(
            "
                #include <netcdf.h>
                int main() {
                int nc_nat_val = NC_NAT;
                return 0;
                }
            "
            HAVE_NETCDF_H)

        target_compile_definitions(cmc PUBLIC WITH_NETCDF=${HAVE_NETCDF_H})
        set(CMC_WITH_NETCDF ON)
        target_include_directories(cmc PUBLIC ${CMC_NETCDF_ROOT_DIR}/include)
        target_link_libraries(cmc ${NETCDF})
        message("-- A netCDF installation has been found")

        if (ENABLE_MPI)
            set(CMAKE_REQUIRED_INCLUDES "${CMC_NETCDF_ROOT_DIR}/include;${MPI_C_INCLUDE_DIRS}")
            set(CMAKE_REQUIRED_LIBRARIES "${NETCDF};MPI::MPI_C")
            set(CMAKE_REQUIRED_LINK_DIRECTORIES "${CMC_NETCDF_ROOT_DIR}/lib")
            check_c_source_compiles(
            "
                #include <netcdf.h>
                #include <netcdf_par.h>
                int main() {
                int nc_independent_val = NC_INDEPENDENT;
                return 0;
                }
            "
            HAVE_NETCDF_PAR_H)

            if (HAVE_NETCDF_PAR_H)
                set(CMC_WITH_NETCDF_PAR ON)    
                target_compile_definitions(cmc PUBLIC WITH_NETCDF_PAR=${HAVE_NETCDF_PAR_H})
                message("-- Found netcdf_par.h")
            else()
                set(CMC_WITH_NETCDF_PAR OFF)  
                message("-- Not Found: netcdf_par.h")
                target_compile_definitions(cmc PUBLIC WITH_NETCDF_PAR=0)
            endif()
        else()
            set(CMC_WITH_NETCDF_PAR OFF)  
            target_compile_definitions(cmc PUBLIC WITH_NETCDF_PAR=0)
        endif()
    else()
        # Find cmake installed netCDF
        find_package(netCDF REQUIRED)
        if(netCDF_FOUND)
            set(CMAKE_REQUIRED_LIBRARIES netCDF::netcdf)
            check_c_source_compiles(
            "
                #include <netcdf.h>
                int main() {
                int nc_nat_val = NC_NAT;
                return 0;
                }
            "
            HAVE_NETCDF_H)
        endif ()
        target_link_libraries(cmc netCDF::netcdf)
        target_compile_definitions(cmc PUBLIC WITH_NETCDF=1)
        set(CMC_WITH_NETCDF ON) 
        message("-- A netCDF installation has been found")

        if (ENABLE_MPI)
            set(CMAKE_REQUIRED_LIBRARIES netCDF::netcdf)
            check_c_source_compiles(
            "
                #include <netcdf.h>
                #include <netcdf_par.h>
                int main() {
                int nc_independent_val = NC_INDEPENDENT;
                return 0;
                }
            "
            HAVE_NETCDF_PAR_H)

            if (HAVE_NETCDF_PAR_H)
                set(CMC_WITH_NETCDF_PAR ON)    
                target_compile_definitions(cmc PUBLIC WITH_NETCDF_PAR=${HAVE_NETCDF_PAR_H})
                message("-- Found netcdf_par.h")
            else()
                set(CMC_WITH_NETCDF_PAR OFF)
                message("-- Not Found: netcdf_par.h") 
                target_compile_definitions(cmc PUBLIC WITH_NETCDF_PAR=0)
            endif()
        else()
            set(CMC_WITH_NETCDF_PAR OFF)  
            target_compile_definitions(cmc PUBLIC WITH_NETCDF_PAR=0)
        endif()
    endif()
else()
    set(CMC_WITH_NETCDF OFF)
    set(CMC_WITH_NETCDF_PAR OFF)
    target_compile_definitions(cmc PUBLIC WITH_NETCDF=0)
    target_compile_definitions(cmc PUBLIC WITH_NETCDF_PAR=0)
endif()

# Check if additional warnings should be activated
if(ENABLE_ALL_WARNINGS)
    target_compile_options(cmc PRIVATE -Wall -Wextra -Wpedantic)
    set(cmc_CXXFLAGS "${cmc_CXXFLAGS} -Wall -Wextra -Wpedantic")
    set(cmc_CFLAGS "${cmc_CFLAGS} -Wall -Wextra -Wpedantic")
    target_compile_definitions(cmc PUBLIC ENABLE_ALL_WARNINGS=1)
    message("-- Compilation flags for all warnings have been set")
    set(CMC_ENABLE_ALL_WARNINGS ON) 
endif()

# Check if t8code will be linked 
if(WITH_T8CODE)
    find_package(T8CODE REQUIRED)
    target_link_libraries(cmc T8CODE::T8)
    message("-- Found t8code installation")
    set(CMC_WITH_T8CODE ON) 
    target_compile_definitions(cmc PUBLIC WITH_T8CODE=1)
else()
    set(CMC_WITH_T8CODE OFF) 
    target_compile_definitions(cmc PUBLIC WITH_T8CODE=0)
endif()

# Define incldue-directories for building and installing
#target_include_directories(cmc PUBLIC "${PROJECT_SOURCE_DIR}")
target_include_directories(cmc PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>
)

# Add this path as an include directory
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Rpath options necessary for shared library install to work correctly in user projects.
set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

# Configure an input file with some given flags
configure_file(CmcConfig.h.in ${PROJECT_BINARY_DIR}/src/cmc_ac_config.h @ONLY)

# Set some compile defintions regarding compialtion flags and versioning
target_compile_definitions(cmc PUBLIC CMC_CC="${CMAKE_C_COMPILER}")
target_compile_definitions(cmc PUBLIC CMC_CFLAGS="${CMC_CFLAGS}")
target_compile_definitions(cmc PUBLIC CMC_CXX="${CMAKE_CXX_COMPILER}")
target_compile_definitions(cmc PUBLIC CMC_CXXFLAGS="${CMC_CXXFLAGS}")
target_compile_definitions(cmc PUBLIC CMC_LDFLAGS="${CMAKE_SHARED_LINKER_FLAGS}")
target_compile_definitions(cmc PUBLIC CMC_PACKAGE_STRING="cmc ${cmc_VERSION}")
target_compile_definitions(cmc PUBLIC CMC_VERSION="${cmc_VERSION}")
target_compile_definitions(cmc PUBLIC CMC_VERSION_MAJOR=${cmc_VERSION_MAJOR})
target_compile_definitions(cmc PUBLIC CMC_VERSION_MINOR=${cmc_VERSION_MINOR})
target_compile_definitions(cmc PUBLIC CMC_VERSION_PATCH=${cmc_VERSION_PATCH})
target_compile_definitions(cmc PUBLIC CMC_VERSION_POINT=${cmc_VERSION_POINT})

# Add the linker flags to the compile definitions
get_target_property(cmc_LIBS_LIST cmc LINK_LIBRARIES)
string(REPLACE ";" " " cmc_LIBS "${cmc_LIBS_LIST}")
target_compile_definitions(cmc PUBLIC cmc_LIBS="${cmc_LIBS}")

# Get in the src-directory
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src)

# Get in the example-directory
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/example)

# Get in the test-directory
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/test)

# Get in the workflows-directory
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/workflows)

# Create config files such that cmc may be found by find_package etc.
include(CMakePackageConfigHelpers)

set(CmcConfigPackageLocation ${CMAKE_INSTALL_LIBDIR}/cmake/cmc)

install(EXPORT CmcTargets
    NAMESPACE cmc::
    DESTINATION ${CmcConfigPackageLocation})

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/CmcConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CmcConfig.cmake
    INSTALL_DESTINATION ${CmcConfigPackageLocation}
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/CmcConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CmcConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/CmcConfigVersion.cmake
        DESTINATION ${CmcConfigPackageLocation})

message("-- Configuration results cmc ${cmc_VERSION_MAJOR}.${cmc_VERSION_MINOR}.${cmc_VERSION_PATCH}:")
message("-- \tBUILD_SHARED_LIBS: \t${BUILD_SHARED_LIBS}")
message("-- \tENABLE_MPI: \t\t${ENABLE_MPI}")
message("-- \tWITH_NETCDF: \t\t${WITH_NETCDF}")
if (WITH_NETCDF)
message("-- \t-- with netcdf_par.h: \t${CMC_WITH_NETCDF_PAR}")
endif()
message("-- \tWITH_T8CODE: \t\t${WITH_T8CODE}")
message("-- \tENABLE_ALL_WARNINGS: \t${ENABLE_ALL_WARNINGS}")
