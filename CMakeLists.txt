###########################################################################
#                 CMake Build File for ICM-Morph
#
#    Written By: Zach Cobell & Eric White
#
###########################################################################
#
# The CMake build system enable ICM-LAVegMod to be deployed and built
# in a cross platform environment. 
#
###########################################################################
#...Set the default build type
IF(DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
        build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
        Release RelWithDebInfo MinSizeRel.")
ELSE()
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
        options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
        RelWithDebInfo MinSizeRel.")
ENDIF()
###########################################################################

###########################################################################
#  GENERAL OPTIONS
###########################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(ICM_LAVegMod Fortran)
ENABLE_LANGUAGE( C )

###########################################################################


###########################################################################
#...Put the static libraries in the CMakeFiles folder so they don't
#   contaminate the build directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles)
###########################################################################


###########################################################################
#...Set the version string (used for libraries)
SET(ICMLVM_VERSION_MAJOR 29)
SET(ICMLVM_VERSION_MINOR  0)
SET(ICMLVM_VERSION_PATCH  0)
SET(ICMLVM_VERSION_STRING ${ICMLVM_VERSION_MAJOR}.${ICMLVM_VERSION_MINOR}.${ICMLVM_VERSION_PATCH})
###########################################################################


###########################################################################
#...Determine architecture specific parameters
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/architecture.cmake)
###########################################################################


###########################################################################
#...Enable/Disable netCDF functionality
OPTION(ENABLE_OUTPUT_NETCDF "Enable netCDF output formats" OFF)
IF(ENABLE_OUTPUT_NETCDF)
    SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
    SET(NETCDF_F90 "YES")
    FIND_PACKAGE(NetCDF)
    IF(NOT NETCDF_FOUND)
        MESSAGE(SEND_ERROR "Could not locate NetCDF. Try setting the command line flag -DNETCDF_DIR=[path]")
    ENDIF()
ENDIF()
###########################################################################


###########################################################################
# ICM_Morph
###########################################################################
SET( LVM LAVegMod_v${ICMLVM_VERSION_STRING} )
ADD_EXECUTABLE( ${LVM}   LVM_acute_salinity.f90
                         LVM_acute_salinity_flt.f90
                         LVM_calc_dispersal_coverage.f90
                         LVM_check_sum.f90                         
                         LVM_land_change.f90
                         LVM_main.f90
                         LVM_mort_est_prob.f90
                         LVM_neighbors.f90
                         LVM_oneway_interp.f90
                         LVM_params.f90
                         LVM_params_alloc.f90
                         LVM_preprocessing.f90
                         LVM_reset_coverages.f90
                         LVM_round_coverages.f90
                         LVM_set_io.f90
                         LVM_sum_unoccupied_flt.f90
                         LVM_sum_unoccupied_land.f90
                         LVM_twoway_interp.f90
                         LVM_update_coverages.f90
                         LVM_update_flotant.f90
                         LVM_weedy_establishment.f90
                         LVM_write_output.f90 )
SET_TARGET_PROPERTIES( ${LVM} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/mod)
SET_TARGET_PROPERTIES( ${LVM} PROPERTIES COMPILE_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG}" )

IF(ENABLE_OUTPUT_NETCDF)
    TARGET_COMPILE_DEFINITIONS( ${LVM} PRIVATE ICM_NETCDF )
    TARGET_INCLUDE_DIRECTORIES( ${LVM} PRIVATE ${NETCDF_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES( ${LVM} PRIVATE ${NETCDF_LIBRARIES})
ENDIF()
