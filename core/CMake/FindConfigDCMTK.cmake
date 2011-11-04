# - find DCMTK libraries and applications
#

##
## THIS FILE WAS MODIFIED TO MATCH THE REQUIREMENTS OF CMTK
##
## - only DCMTK libraries are tested and configured that are actually used by CMTK
##
##

#  DCMTK_INCLUDE_DIRS   - Directories to include to use DCMTK
#  DCMTK_LIBRARIES     - Files to link against to use DCMTK
#  DCMTK_FOUND         - If false, don't try to use DCMTK
#  DCMTK_DIR           - (optional) Source directory for DCMTK
#  DCMTK_DCMDICTPATH   - Path to dicom.dic data dictionary
#
# DCMTK_DIR can be used to make it simpler to find the various include
# directories and compiled libraries if you've just compiled it in the
# source tree. Just set it to the root of the tree where you extracted
# the source (default to /usr/include/dcmtk/)

#=============================================================================
# Copyright 2004-2009 Kitware, Inc.
# Copyright 2009-2010 Mathieu Malaterre <mathieu.malaterre@gmail.com>
# Copyright 2010 Thomas Sondergaard <ts@medical-insight.com>
# Copyright 2011 SRI International
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

#
# Written for VXL by Amitha Perera.
# Upgraded for GDCM by Mathieu Malaterre.
# Modified for EasyViz by Thomas Sondergaard.
# 
# Wed Jul 27 2011 Ankur Sinha <ankursinha AT fedoraproject DOT org> 
# - Add all dcmtk libs
# - Add usr/lib to paths
#

if(NOT DCMTK_FOUND AND NOT DCMTK_DIR)
  set(DCMTK_DIR
    "/usr/include/dcmtk/"
    CACHE
    PATH
    "Root of DCMTK source tree (optional).")
  mark_as_advanced(DCMTK_DIR)
endif()


foreach(lib
    dcmdata
##    dcmdsig
    dcmimage
    dcmimgle
    dcmjpeg
    dcmjpls
##    dcmnet
##    dcmpstat
##    dcmqrdb
##    dcmsr
##    dcmtls
##    dcmwlm
    ijg12
    ijg16
    ijg8
##    libi2d
    oflog
    ofstd)



  find_library(DCMTK_${lib}_LIBRARY
    ${lib}
    PATHS
    ${DCMTK_DIR}/${lib}/libsrc
    ${DCMTK_DIR}/${lib}/libsrc/Release
    ${DCMTK_DIR}/${lib}/libsrc/Debug
    ${DCMTK_DIR}/${lib}/Release
    ${DCMTK_DIR}/${lib}/Debug
    ${DCMTK_DIR}/lib
    /usr/lib/dcmtk)

  mark_as_advanced(DCMTK_${lib}_LIBRARY)

  if(DCMTK_${lib}_LIBRARY)
    list(APPEND DCMTK_LIBRARIES ${DCMTK_${lib}_LIBRARY})
  endif()

endforeach()


set(DCMTK_config_TEST_HEADER osconfig.h)
set(DCMTK_dcmdata_TEST_HEADER dctypes.h)
set(DCMTK_dcmimage_TEST_HEADER dicoimg.h)
set(DCMTK_dcmimgle_TEST_HEADER dcmimage.h)
set(DCMTK_dcmjpeg_TEST_HEADER djdecode.h)
set(DCMTK_dcmjpls_TEST_HEADER djcodecd.h)
##set(DCMTK_dcmnet_TEST_HEADER assoc.h)
##set(DCMTK_dcmpstat_TEST_HEADER dcmpstat.h)
##set(DCMTK_dcmqrdb_TEST_HEADER dcmqrdba.h)
##set(DCMTK_dcmsign_TEST_HEADER sicert.h)
##set(DCMTK_dcmsr_TEST_HEADER dsrtree.h)
##set(DCMTK_dcmtls_TEST_HEADER tlslayer.h)
##set(DCMTK_dcmwlm_TEST_HEADER wldsfs.h)
set(DCMTK_ofstd_TEST_HEADER ofstdinc.h)
set(DCMTK_oflog_TEST_HEADER oflog.h)

foreach(dir
    config
    dcmdata
    dcmimage
    dcmimgle
    dcmjpeg
    dcmjpls
##    dcmnet
##    dcmpstat
##    dcmqrdb
##    dcmsign
##    dcmsr
##    dcmtls
##    dcmwlm
    ofstd
    oflog)

  find_path(DCMTK_${dir}_INCLUDE_DIR
    ${DCMTK_${dir}_TEST_HEADER}
    PATHS
    ${DCMTK_DIR}/${dir}/include
    ${DCMTK_DIR}/${dir}
    ${DCMTK_DIR}/include/${dir}
    /usr/include/dcmtk)

  mark_as_advanced(DCMTK_${dir}_INCLUDE_DIR)

  if(DCMTK_${dir}_INCLUDE_DIR)
    list(APPEND
      DCMTK_INCLUDE_DIRS
      ${DCMTK_${dir}_INCLUDE_DIR})
  endif()
endforeach()

if(WIN32)
  list(APPEND DCMTK_LIBRARIES netapi32 wsock32)
endif()

if(DCMTK_ofstd_INCLUDE_DIR)
  get_filename_component(DCMTK_dcmtk_INCLUDE_DIR
    ${DCMTK_ofstd_INCLUDE_DIR}
    PATH
    CACHE)
  list(APPEND DCMTK_INCLUDE_DIRS ${DCMTK_dcmtk_INCLUDE_DIR})
  mark_as_advanced(DCMTK_dcmtk_INCLUDE_DIR)
endif()

include(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)
find_package_handle_standard_args(DCMTK DEFAULT_MSG
  DCMTK_config_INCLUDE_DIR
  DCMTK_ofstd_INCLUDE_DIR
  DCMTK_ofstd_LIBRARY
  DCMTK_dcmdata_INCLUDE_DIR
  DCMTK_dcmdata_LIBRARY
  DCMTK_dcmimgle_INCLUDE_DIR
  DCMTK_dcmimgle_LIBRARY)

# Compatibility: This variable is deprecated
set(DCMTK_INCLUDE_DIR ${DCMTK_INCLUDE_DIRS})

##foreach(executable dcmdump dcmdjpeg dcmdrle)
##  string(TOUPPER ${executable} EXECUTABLE)
##  find_program(DCMTK_${EXECUTABLE}_EXECUTABLE ${executable} ${DCMTK_DIR}/bin)
##  mark_as_advanced(DCMTK_${EXECUTABLE}_EXECUTABLE)
##endforeach()

FIND_PATH(DCMTK_DCMDICTPATH
  dicom.dic
  PATHS
  ${DCMTK_DIR}/lib
  ${DCMTK_DIR}/share
  /usr/include/dcmtk
  /usr/share/dcmtk)

# The following using pieces contributed by Yaroslav Halchenko based on CMake's own CHECK_CXX_SOURCE_COMPILES
FUNCTION(CheckLibraryDependency _required _lib _key)
  IF(NOT DEFINED REQUIRED_LIBRARY_${_lib})
    SET(CHECK_CXX_SOURCE_COMPILES_ADD_LIBRARIES "-DLINK_LIBRARIES:STRING=${${_required}}")
    FILE(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx" "int main(void) { return 0; }\n")
    
    #MESSAGE( "D: Checking for ${_lib}" )
    TRY_COMPILE(_VAR_IGNORE_
      ${CMAKE_BINARY_DIR}
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx
      COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
      CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
      "${CHECK_CXX_SOURCE_COMPILES_ADD_LIBRARIES}"
      OUTPUT_VARIABLE OUTPUT)
    
    #MESSAGE( "D: Output was ${OUTPUT}" )
    IF("${OUTPUT}" MATCHES ".*undefined reference to.*${_key}.*")
      MESSAGE(STATUS "+ required library ${_lib}" )
      SET(REQUIRED_LIBRARY_${_lib} TRUE CACHE BOOL "Library ${_lib} is required to resolve dependencies.")
    ELSE()
      MESSAGE(STATUS "- tested library ${_lib}" )
      SET(REQUIRED_LIBRARY_${_lib} FALSE CACHE BOOL "Library ${_lib} is not required to resolve dependencies.")
    ENDIF()
  ENDIF(NOT DEFINED REQUIRED_LIBRARY_${_lib})
  
  IF( REQUIRED_LIBRARY_${_lib} )  
    SET(${_required} ${${_required}};${_lib} PARENT_SCOPE)
  ENDIF( REQUIRED_LIBRARY_${_lib} )  
ENDFUNCTION()

    # Detect missing DCMTK library dependencies by testing the "usual suspects"
##    CheckLibraryDependency(DCMTK_LIBRARIES wrap hosts_access)
CheckLibraryDependency(DCMTK_LIBRARIES png png_write_image)
CheckLibraryDependency(DCMTK_LIBRARIES tiff TIFFGetVersion)
CheckLibraryDependency(DCMTK_LIBRARIES CharLS JpegLsReadHeader)    
##    CheckLibraryDependency(DCMTK_LIBRARIES xml2 xmlGetProp)

