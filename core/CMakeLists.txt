##
##  Copyright 1997-2010 Torsten Rohlfing
##
##  Copyright 2004-2011 SRI International
##
##  This file is part of the Computational Morphometry Toolkit.
##
##  http://www.nitrc.org/projects/cmtk/
##
##  The Computational Morphometry Toolkit is free software: you can
##  redistribute it and/or modify it under the terms of the GNU General Public
##  License as published by the Free Software Foundation, either version 3 of
##  the License, or (at your option) any later version.
##
##  The Computational Morphometry Toolkit is distributed in the hope that it
##  will be useful, but WITHOUT ANY WARRANTY; without even the implied
##  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with the Computational Morphometry Toolkit.  If not, see
##  <http://www.gnu.org/licenses/>.
##
##  $Revision$
##
##  $LastChangedDate$
##
##  $LastChangedBy$
##

##
## This file borrows heavily from the analogous InsightToolkit file
##

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW) ## libraries with full paths vs. search paths
ENDIF(COMMAND cmake_policy)
SET(CMAKE_OVERRIDE_COMPILER_MISMATCH 1)

PROJECT(CMTK)
SET(CMTK_VERSION_MAJOR "2")
SET(CMTK_VERSION_MINOR "1")
SET(CMTK_VERSION_PATCH "0svn")

SET(CMTK_VERSION_STRING "${CMTK_VERSION_MAJOR}.${CMTK_VERSION_MINOR}.${CMTK_VERSION_PATCH}")

#-----------------------------------------------------------------------------
# Test for some required system information.
INCLUDE (${CMAKE_ROOT}/Modules/CMakeDetermineSystem.cmake)
INCLUDE (${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeSystem.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityC.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityCXX.cmake)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake)

#-----------------------------------------------------------------------------
# Output directories.
IF(NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${CMTK_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
ENDIF(NOT LIBRARY_OUTPUT_PATH)
IF(NOT EXECUTABLE_OUTPUT_PATH)
  SET (EXECUTABLE_OUTPUT_PATH ${CMTK_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
ENDIF(NOT EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
IF(NOT CMTK_BUILD_LIB_DIR) # variable used for generating CMTKConfig.cmake for the build tree
  SET(CMTK_BUILD_LIB_DIR "/bin")
ENDIF(NOT CMTK_BUILD_LIB_DIR)

SET(CMTK_LIBRARY_PATH "${LIBRARY_OUTPUT_PATH}")
SET(CMTK_EXECUTABLE_PATH "${EXECUTABLE_OUTPUT_PATH}")
SET(CXX_TEST_PATH ${EXECUTABLE_OUTPUT_PATH})

# set CMTK_DIR so it can be used by subprojects
SET(CMTK_DIR "${CMAKE_BINARY_DIR}" CACHE INTERNAL "CMTK dir to be used by subprojects")

#-----------------------------------------------------------------------------
# Setup install target directories
IF(NOT CMTK_INSTALL_BIN_DIR)
  SET(CMTK_INSTALL_BIN_DIR bin CACHE PATH "Binary files will be installed here")
ENDIF(NOT CMTK_INSTALL_BIN_DIR)

IF(NOT CMTK_INSTALL_LIB_DIR)
  # 64 bit libraries are not installed in the same directories on solaris
  # Ref: http://docs.sun.com/app/docs/doc/816-5138/6mba6ua5m?a=view
  SET(sun64 )
  IF(CMAKE_SYSTEM MATCHES "SunOS.*" AND CMAKE_SIZEOF_VOID_P EQUAL 8)
    EXEC_PROGRAM("uname"
      ARGS "-p"
      OUTPUT_VARIABLE sunproc
    )
    IF("${sunproc}" STREQUAL "sparc")
      SET(sun64 "/sparcv9")
    ENDIF("${sunproc}" STREQUAL "sparc")

    IF("${sunproc}" STREQUAL "i386")
      SET(sun64 "/amd64")
    ENDIF("${sunproc}" STREQUAL "i386")

    IF("${sun64}" STREQUAL "")
      # sun64 is not set - solaris has been ported to a new architecture?
      MESSAGE("Warning: Unknown processor '${sunproc}'. The libraries may not be installed in the right place.")
    ENDIF("${sun64}" STREQUAL "")
  ENDIF(CMAKE_SYSTEM MATCHES "SunOS.*" AND CMAKE_SIZEOF_VOID_P EQUAL 8)

  SET(CMTK_INSTALL_LIB_DIR "lib/cmtk${sun64}" CACHE PATH "Library files will be installed here")
ENDIF(NOT CMTK_INSTALL_LIB_DIR)

IF(NOT CMTK_INSTALL_DATA_DIR)
  SET(CMTK_INSTALL_DATA_DIR "share/cmtk" CACHE PATH "Data files will be installed here")
ENDIF(NOT CMTK_INSTALL_DATA_DIR)

IF(NOT CMTK_INSTALL_INCLUDE_DIR)
  SET(CMTK_INSTALL_INCLUDE_DIR "include/cmtk" CACHE PATH "Data files will be installed here")
ENDIF(NOT CMTK_INSTALL_INCLUDE_DIR)

#-----------------------------------------------------------------------------
# Check for system headers.
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES (dirent.h HAVE_DIRENT_H)
CHECK_INCLUDE_FILES (execinfo.h HAVE_EXECINFO_H)
CHECK_INCLUDE_FILES (fcntl.h HAVE_FCNTL_H)
CHECK_INCLUDE_FILES (ieeefp.h HAVE_IEEEFP_H)
CHECK_INCLUDE_FILES (inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES (malloc.h HAVE_MALLOC_H)
CHECK_INCLUDE_FILES (pthread.h HAVE_PTHREAD_H)
CHECK_INCLUDE_FILES (stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILES (termios.h HAVE_TERMIOS_H)
CHECK_INCLUDE_FILES (unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILES (values.h HAVE_VALUES_H)

CHECK_INCLUDE_FILES (sys/ioctl.h HAVE_SYS_IOCTL_H)
CHECK_INCLUDE_FILES (sys/procfs.h HAVE_SYS_PROCFS_H)
CHECK_INCLUDE_FILES (sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES (sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES (sys/times.h HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILES (sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES (sys/utsname.h HAVE_SYS_UTSNAME_H)

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(WORDS_BIGENDIAN)

#-----------------------------------------------------------------------------
# Test NaN high bit
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMake/TestQnanhibit.cmake)
TEST_QNANHIBIT(QNANHIBIT_VALUE ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
IF(QNANHIBIT_VALUE)
  SET(QNANHIBIT 1 CACHE INTERNAL "The 22nd bit of 32-bit floating-point quiet NaN.")
ELSE(QNANHIBIT_VALUE)
  SET(QNANHIBIT 0 CACHE INTERNAL "The 22nd bit of 32-bit floating-point quiet NaN.")
ENDIF(QNANHIBIT_VALUE)

#-----------------------------------------------------------------------------
# Check for at least one of the typical STL hash map implementations

INCLUDE (CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(unordered_map HAVE_UNORDERED_MAP)
IF(NOT HAVE_UNORDERED_MAP)
  CHECK_INCLUDE_FILE_CXX(tr1/unordered_map HAVE_UNORDERED_MAP_TR1)
  IF(NOT HAVE_UNORDERED_MAP_TR1)
    CHECK_INCLUDE_FILE_CXX(hash_map HAVE_HASH_MAP)
    IF(NOT HAVE_HASH_MAP)
      CHECK_INCLUDE_FILE_CXX(hash_map.h HAVE_HASH_MAP_H)
    ENDIF(NOT HAVE_HASH_MAP)
  ENDIF(NOT HAVE_UNORDERED_MAP_TR1)
ENDIF(NOT HAVE_UNORDERED_MAP)

# set a single flag is we have ANY hash map implementation in STL
SET(HAVE_STL_HASH_MAP ${HAVE_UNORDERED_MAP}${HAVE_UNORDERED_MAP_TR1}${HAVE_HASH_MAP}${HAVE_HASH_MAP_H} CACHE INTERNAL "Flag for presence of ANY hash map implementation in the C++ STL")

#-----------------------------------------------------------------------------
# build configuration options.
OPTION(BUILD_SHARED_LIBS "Build CMTK with shared libraries." OFF)
SET(CMTK_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

OPTION(CMTK_BUILD_UNSTABLE "Build 'unstable' library and tools" OFF)
MARK_AS_ADVANCED(CMTK_BUILD_UNSTABLE)

OPTION(CMTK_BUILD_STACKTRACE "Build tools with stack trace printing in case of crash (increases build time)" ON)
MARK_AS_ADVANCED(CMTK_BUILD_STACKTRACE)

OPTION(CMTK_BUILD_DEMO "Build demonstration code. Strange things will happen if you turn this on." OFF)
MARK_AS_ADVANCED(CMTK_BUILD_DEMO)
IF(CMTK_BUILD_DEMO)
  MESSAGE( WARNING "Demo mode is enabled. Make sure not to use this build for production use!")
ENDIF(CMTK_BUILD_DEMO)

#-----------------------------------------------------------------------------
# Add an option to use or not to use MPI distributed-memory parallelism

OPTION(CMTK_USE_MPI "Use MPI for distributed-memory parallelism" OFF)
MARK_AS_ADVANCED(CMTK_USE_MPI)

OPTION(CMTK_USE_SMP "Enable shared-memory parallelism (requires at least POSIX Threads, OpenMP, or Windows threads)" ON)

#-----------------------------------------------------------------------------
# Add an option to enable or disable POSIX threads if pthread.h header exists

IF(HAVE_PTHREAD_H)
  OPTION(CMTK_USE_PTHREADS "Use POSIX Threads for SMP parallelism" ON)
  MARK_AS_ADVANCED(CMTK_USE_PTHREADS)

  IF(NOT CMTK_USE_SMP)
    SET(CMTK_USE_PTHREADS OFF)
  ENDIF(NOT CMTK_USE_SMP)
ENDIF(HAVE_PTHREAD_H)

#-----------------------------------------------------------------------------
# for the gnu compiler a -D_PTHREADS is needed on sun
# for the native compiler a -mt flag is needed on the sun
IF(CMTK_USE_PTHREADS)

  IF(${CMAKE_SYSTEM} MATCHES "SunOS.*")
    IF(CMAKE_COMPILER_IS_GNUCXX)
      SET(CMTK_REQUIRED_CXX_FLAGS "${CMTK_REQUIRED_CXX_FLAGS} -D_PTHREADS")
    ELSE(CMAKE_COMPILER_IS_GNUCXX)
      SET(CMTK_REQUIRED_CXX_FLAGS "${CMTK_REQUIRED_CXX_FLAGS} -mt")
      SET(CMTK_REQUIRED_C_FLAGS "${CMTK_REQUIRED_C_FLAGS} -mt")
    ENDIF(CMAKE_COMPILER_IS_GNUCXX)
  ENDIF(${CMAKE_SYSTEM} MATCHES "SunOS.*")
ENDIF(CMTK_USE_PTHREADS)

# mingw thread support
IF(MINGW)
  SET(CMTK_REQUIRED_CXX_FLAGS "${CMTK_REQUIRED_CXX_FLAGS} -mthreads")
  SET(CMTK_REQUIRED_C_FLAGS "${CMTK_REQUIRED_C_FLAGS} -mthreads")
  SET(CMTK_REQUIRED_LINK_FLAGS "${CMTK_REQUIRED_LINK_FLAGS} -mthreads")
ENDIF(MINGW)

#-----------------------------------------------------------------------------
# Check for, and configure, OpenMP

IF(NOT DEFINED CMTK_USE_OPENMP)
  IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION} STRLESS 2.6.2)
    MESSAGE("Warning.  We need at least CMake 2.6.2 to check for OpenMP support. You have version ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}")
  ELSE (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION} STRLESS 2.6.2)
    FIND_PACKAGE(OpenMP)
    IF(OPENMP_FOUND)
      OPTION(CMTK_USE_OPENMP "Use OpenMP for SMP parallelism" ON)
      MARK_AS_ADVANCED(CMTK_USE_OPENMP)
      
      IF(NOT CMTK_USE_SMP)
	SET(CMTK_USE_OPENMP OFF)
      ENDIF(NOT CMTK_USE_SMP)
    ENDIF(OPENMP_FOUND)
  ENDIF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION} STRLESS 2.6.2)
ENDIF(NOT DEFINED CMTK_USE_OPENMP)

IF(CMTK_USE_OPENMP)
  SET(CMTK_REQUIRED_CXX_FLAGS "${CMTK_REQUIRED_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  SET(CMTK_REQUIRED_LINK_FLAGS "${CMTK_REQUIRED_LINK_FLAGS} ${OpenMP_C_FLAGS}")
ENDIF(CMTK_USE_OPENMP)

#-----------------------------------------------------------------------------
# Check for, and configure, Grand Central Dispatch

CHECK_INCLUDE_FILES(dispatch/dispatch.h HAVE_DISPATCH_H)
IF(HAVE_DISPATCH_H)
  OPTION(CMTK_USE_GCD "Use Grand Central Dispatch for SMP parallelism with system-level load balancing" ON)
  
  IF(NOT CMTK_USE_SMP)
    SET(CMTK_USE_GCD OFF)
  ENDIF(NOT CMTK_USE_SMP)
  
  IF(CMTK_USE_GCD)
  ENDIF(CMTK_USE_GCD)
ENDIF(HAVE_DISPATCH_H)

#-----------------------------------------------------------------------------
# Check for, and configure, OpenCL

##    OPTION(CMTK_USE_OPENCL "Use OpenCL for GPU computing" OFF)
##  IF(CMTK_USE_OPENCL)
##  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMake/FindOpenCL.cmake)
##  IF(OPENCL_FOUND)
##  ENDIF(OPENCL_FOUND)
##
##    ADD_DEFINITIONS(-DCMTK_USE_OPENCL)
##  ENDIF(CMTK_USE_OPENCL)

#-----------------------------------------------------------------------------
# Check for, and configure, CUDA

OPTION(CMTK_USE_CUDA "Use CUDA for GPU computing" OFF)

IF(CMTK_USE_CUDA)

  # check using cmake's own FindCUDA above v2.8 or use our own for earlier versions
  IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} STRLESS 2.8)
    INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMake/FindCUDA.cmake)
  ELSE (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} STRLESS 2.8)
    FIND_PACKAGE(CUDA)
  ENDIF(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} STRLESS 2.8)
  
ENDIF(CMTK_USE_CUDA)

#-----------------------------------------------------------------------------
# Check whether compiler supports variable-size automatic arrays, i.e.,
# things like "const int size = 10; const char string[size];" rather than
# "const char string[10];"

INCLUDE(${CMAKE_ROOT}/Modules/CheckCXXSourceCompiles.cmake)
CHECK_CXX_SOURCE_COMPILES("class A {}; class B : public A {}; int main() { int size = 10; B array[size]; }" CMTK_COMPILER_VAR_AUTO_ARRAYSIZE)

#-----------------------------------------------------------------------------
# Configure Dart testing support.  This should be done before any
# MESSAGE(FATAL_ERROR ...) commands are invoked.

INCLUDE(${CMAKE_ROOT}/Modules/CTest.cmake)
MARK_AS_ADVANCED(TCL_TCLSH DART_ROOT)
ENABLE_TESTING()
IF(BUILD_TESTING)
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)

  OPTION(CMTK_TESTING_MEMORYCHECK "Test using memory checking" OFF)
  MARK_AS_ADVANCED(CMTK_TESTING_MEMORYCHECK)

  IF(NOT CMTK_DATA_ROOT)
    GET_FILENAME_COMPONENT(CMTK_SVN_ROOT ${CMAKE_CURRENT_SOURCE_DIR} PATH)
    SET(CMTK_DATA_ROOT ${CMTK_SVN_ROOT}/data CACHE PATH "Root directory of the CMTK data tree")
  ENDIF(NOT CMTK_DATA_ROOT)
  MARK_AS_ADVANCED(CMTK_DATA_ROOT)  
ENDIF(BUILD_TESTING)

SET(CMTK_ROOT_PATH_SRI24 "CMTK_ROOT_PATH_SRI24-NOTFOUND" CACHE PATH "Root path of the SRI24 atlas data (usually ends with \"sri24\")" )

#-----------------------------------------------------------------------------
# Configure Doxygen support

FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
  OPTION(BUILD_DOCUMENTATION "Build html-based library API documentation using Doxygen" OFF)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  IF(BUILD_DOCUMENTATION)
    ADD_CUSTOM_TARGET(doc COMMAND doxygen DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ DESTINATION ${CMTK_INSTALL_DATA_DIR}/doc/html/ COMPONENT documentation)
  ENDIF(BUILD_DOCUMENTATION)
ENDIF(DOXYGEN_FOUND)

#-----------------------------------------------------------------------------
# Add an option to use double-precision floats for coordinates.

OPTION(CMTK_COORDINATES_DOUBLE "Use double precision floats for coordinates (turn off to save memory)" ON)
MARK_AS_ADVANCED(CMTK_COORDINATES_DOUBLE) 

OPTION(CMTK_DATA_DOUBLE "Use double precision floats for data exchange between typed arrays" ON)
MARK_AS_ADVANCED(CMTK_DATA_DOUBLE) 

OPTION(CMTK_NUMERICS_DOUBLE "Use double precision floats for numerical algorithms" ON)
MARK_AS_ADVANCED(CMTK_NUMERICS_DOUBLE)

#-----------------------------------------------------------------------------
# Add an option to use or not use zlib library

OPTION(CMTK_BUILD_ZLIB "Build bundled zlib library, even if system zlib exists" OFF)
FIND_PACKAGE(ZLIB)
IF(CMTK_BUILD_ZLIB)
  MESSAGE( STATUS "Building bundled zlib by user configuration" )
  SET(CMTK_BUILD_ZLIB 1)
  
  SET(ZLIB_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Utilities/zlib")
  SET(ZCONF_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/Utilities/zlib")
  SET(ZLIB_LIBRARIES "cmtkZlib" CACHE INTERNAL "")
  INCLUDE_DIRECTORIES(BEFORE ${ZCONF_INCLUDE_DIR}  ${ZLIB_INCLUDE_DIR})
ELSE(CMTK_BUILD_ZLIB)
  IF(ZLIB_FOUND)
    MESSAGE( STATUS "Using system zlib" )
  ELSE(ZLIB_FOUND)
    MESSAGE( STATUS "Building bundled zlib due to missing system zlib" )
    SET(CMTK_BUILD_ZLIB 1)
    
    SET(ZLIB_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Utilities/zlib")
    SET(ZCONF_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/Utilities/zlib")
    SET(ZLIB_LIBRARIES "cmtkZlib" CACHE INTERNAL "")
    INCLUDE_DIRECTORIES(BEFORE ${ZCONF_INCLUDE_DIR} ${ZLIB_INCLUDE_DIR})
  ENDIF(ZLIB_FOUND)
ENDIF(CMTK_BUILD_ZLIB)

#-----------------------------------------------------------------------------
# Add an option to use or not use mxml library

OPTION(CMTK_BUILD_MXML "Build bundled mxml library, even if system library exists" OFF)

FIND_PACKAGE(MXML)
IF(CMTK_BUILD_MXML)
  MESSAGE( STATUS "Building bundled mxml by user configuration" )
  SET(CMTK_BUILD_MXML 1)
  
  SET(MXML_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Utilities/mxml")
  INCLUDE_DIRECTORIES(BEFORE ${MXML_INCLUDE_DIR})
  SET(MXML_LIBRARIES "cmtkMxml" CACHE INTERNAL "")
ELSE(CMTK_BUILD_MXML)
  IF(MXML_FOUND)
    MESSAGE( STATUS "Using system mxml" )
  ELSE(MXML_FOUND)
    MESSAGE( STATUS "Building bundled mxml due to missing system mxml" )
    SET(CMTK_BUILD_MXML 1)
    
    SET(MXML_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Utilities/mxml")
    INCLUDE_DIRECTORIES(BEFORE ${MXML_INCLUDE_DIR})
    SET(MXML_LIBRARIES "cmtkMxml" CACHE INTERNAL "")
  ENDIF(MXML_FOUND)
ENDIF(CMTK_BUILD_MXML)

#-----------------------------------------------------------------------------
# Add an option to use or not use bzip2 library

FIND_PACKAGE(BZip2)
IF(BZIP2_FOUND)
  MESSAGE( STATUS "Using system libbz2" )
  SET(CMTK_BZIP2_LIBS ${BZIP2_LIBRARIES})
  OPTION(CMTK_USE_BZIP2 "Compile support for on-the-fly BZip2 decompression" ON)
  MARK_AS_ADVANCED(CMTK_USE_BZIP2)
ENDIF(BZIP2_FOUND)

#-----------------------------------------------------------------------------
# Add an option to use or not use lzma decompression library

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMake/FindLZMA.cmake)
IF(LZMA_FOUND)
  MESSAGE( STATUS "Using system lzma library" )
  SET(CMTK_LZMA_LIBS ${LZMA_LIBRARY})
  OPTION(CMTK_USE_LZMA "Compile support for on-the-fly LZMA decompression" ON)
  MARK_AS_ADVANCED(CMTK_USE_LZMA)
ENDIF(LZMA_FOUND)

IF(CMTK_USE_LZMA)
  INCLUDE_DIRECTORIES(${LZMA_INCLUDE_DIR})
ENDIF(CMTK_USE_LZMA)

#-----------------------------------------------------------------------------
# Add an option to build with DICOM support

OPTION(CMTK_USE_DCMTK "Build with optional DICOM support" ON)

IF(CMTK_USE_DCMTK)
  FIND_PACKAGE(DCMTK)

  OPTION(CMTK_BUILD_DCMTK "Build bundled DCMTK library, even if system library exists" OFF)

  IF(NOT DCMTK_FOUND)
    MESSAGE( STATUS "DCMTK library not found - defaulting to bundled version." )
    SET(CMTK_BUILD_DCMTK ON CACHE BOOL "Build bundled DCMTK library, even if system library exists" FORCE )
  ENDIF(NOT DCMTK_FOUND)

  IF(CMTK_BUILD_DCMTK)
    MESSAGE( STATUS "Building bundled DCMTK library for DICOM support." )
    SET(DCMTK_SOURCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/Utilities/dcmtk)
    SET(DCMTK_BINARY_ROOT ${CMAKE_CURRENT_BINARY_DIR}/Utilities/dcmtk)
    SET(DCMTK_INCLUDE_DIR "${DCMTK_SOURCE_ROOT}")
    SET(DCMTK_INCLUDE_DIRS "${DCMTK_BINARY_ROOT};${DCMTK_INCLUDE_DIR};${DCMTK_INCLUDE_DIR}/dcmtk/dcmjpeg;${DCMTK_INCLUDE_DIR}/dcmtk/dcmimage;${DCMTK_INCLUDE_DIR}/dcmtk/dcmimgle;${DCMTK_INCLUDE_DIR}/dcmtk/dcmdata;${DCMTK_INCLUDE_DIR}/dcmtk/ofstd")
    INCLUDE_DIRECTORIES(BEFORE ${DCMTK_INCLUDE_DIRS})
    SET(DCMTK_LIBRARIES "D4CMTK" CACHE INTERNAL "The bundled DCMTK library.")
    SET(CMAKE_DCMDICTPATH ${DCMTK_SOURCE_ROOT})
    
    SET(DCM_DICT_DEFAULT_PATH "${CMAKE_DCMDICTPATH}")
  ELSE(CMTK_BUILD_DCMTK)

    MESSAGE( STATUS "Using system DCMTK library." )

    # DCMTK needs png and tiff libraries but FindDCMTK doesn't provide dependencies
    SET(DCMTK_LIBRARIES_EXTRA png tiff CharLS xml2 ssl)
    
  ENDIF(CMTK_BUILD_DCMTK)
  
  IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    ADD_DEFINITIONS(-DHAVE_CONFIG_H)
    SET(CMTK_REQUIRED_CXX_FLAGS "${CMTK_REQUIRED_CXX_FLAGS} -DHAVE_CONFIG_H")
  ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    
ELSE(CMTK_USE_DCMTK)

  MESSAGE( STATUS "Building without DICOM support." )
  SET(CMTK_BUILD_DCMTK OFF)
  UNSET(DCMTK_INCLUDE_DIRS CACHE)
  UNSET(DCMTK_LIBRARIES CACHE)

ENDIF(CMTK_USE_DCMTK)

#-----------------------------------------------------------------------------
# Add an option to use or not use sqlite3 library

OPTION(CMTK_USE_SQLITE "Use sqlite3 library for optional database support" OFF)
IF(CMTK_USE_SQLITE)
  OPTION(CMTK_BUILD_SQLITE "Build bundled sqlite3 library, even if system library exists" OFF)
  IF(CMTK_BUILD_SQLITE)
    MESSAGE( STATUS "Building bundled sqlite3 library by user configuration" )
    SET(CMTK_BUILD_SQLITE 1)
    
    SET(SQLITE_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Utilities/sqlite3")
    SET(CMTK_SQLITE_LIB "cmtksqlite3" CACHE INTERNAL "")
    INCLUDE_DIRECTORIES(BEFORE ${SQLITE_INCLUDE_DIR})
  ELSE(CMTK_BUILD_SQLITE)
    INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMake/FindSQLite3.cmake)
    IF(SQLITE3_FOUND AND SQLITE3_VERSION STRGREATER "3.5.0")
      MESSAGE( STATUS "Using system sqlite3 library" )
      SET(SQLITE_INCLUDE_DIR "${SQLITE3_INCLUDE_DIR}")
      SET(CMTK_SQLITE_LIB ${SQLITE3_LIBRARIES} CACHE INTERNAL "")
    ELSE(SQLITE3_FOUND AND SQLITE3_VERSION STRGREATER "3.5.0")
      MESSAGE( STATUS "Building bundled sqlite3 library due to missing or outdated system library" )
      SET(CMTK_BUILD_SQLITE 1)
      
      SET(SQLITE_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Utilities/sqlite3")
      SET(CMTK_SQLITE_LIB "cmtksqlite3" CACHE INTERNAL "")
      INCLUDE_DIRECTORIES(BEFORE ${SQLITE_INCLUDE_DIR})
    ENDIF(SQLITE3_FOUND AND SQLITE3_VERSION STRGREATER "3.5.0")
  ENDIF(CMTK_BUILD_SQLITE)
ENDIF(CMTK_USE_SQLITE)

#-----------------------------------------------------------------------------
# Add an option to use or not use TEEM (for NrrdIO support)

OPTION(CMTK_BUILD_NRRD "Build NrrdIO library for NRRD import/export" ON)

IF(CMTK_BUILD_NRRD)
  ADD_DEFINITIONS(-DTEEM_ZLIB)
  SET(CMTK_REQUIRED_CXX_FLAGS "${CMTK_REQUIRED_CXX_FLAGS} -DTEEM_ZLIB")
  SET(NRRD_LIBRARIES "NrrdIO" CACHE INTERNAL "The library name for the NrrdIO library")
  SET(NRRD_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/Utilities/NrrdIO ${CMAKE_CURRENT_SOURCE_DIR}/Utilities/NrrdIO)
ENDIF(CMTK_BUILD_NRRD)

#-----------------------------------------------------------------------------
# Set up CPack support, to be able to distribute binary packages.
#

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cmtkConfigureCPack.cmake)

#-----------------------------------------------------------------------------
# Configure bundled utilities.

ADD_SUBDIRECTORY(Utilities)

#-----------------------------------------------------------------------------
# Add an option to enable or disable QT support and GUI applications

IF(NOT CMTK_USE_MPI)
  OPTION(CMTK_USE_QT "Use Qt (GUI applications need this)" OFF)
  IF(CMTK_USE_QT)
    SET (QT_MT_REQUIRED true)
    SET (QT_MIN_VERSION "4.3.0")
    FIND_PACKAGE(Qt4 COMPONENTS QtCore QtGui REQUIRED)
    IF (QT_FOUND)
      MESSAGE("Qt4 Found OK (${qt_version_str})")
      INCLUDE(${QT_USE_FILE})
    ELSE(QT_FOUND)
      MESSAGE(FATAL_ERROR "No Qt4")
    ENDIF(QT_FOUND)
  ENDIF(CMTK_USE_QT)
ENDIF(NOT CMTK_USE_MPI)

#-----------------------------------------------------------------------------
# Add options to enable or disable different compoments

OPTION(BUILD_APPS "Build command line tools" ON)
IF(BUILD_APPS)
  SUBDIRS(apps)
ENDIF(BUILD_APPS)

IF(NOT CMTK_USE_MPI)
  IF(CMTK_USE_QT)
    OPTION(BUILD_GUI "Build GUI applications" ON)
  ENDIF(CMTK_USE_QT)
ENDIF(NOT CMTK_USE_MPI)

OPTION(BUILD_TESTING "Build test code" ON)
IF (BUILD_TESTING)
  SUBDIRS(testing)
ENDIF (BUILD_TESTING)

#-----------------------------------------------------------------------------
# Create the list of include directories needed for header files.
INCLUDE(${CMTK_SOURCE_DIR}/cmtkIncludeDirectories.cmake)

INCLUDE_DIRECTORIES(
  ${CMTK_INCLUDE_DIRS_BUILD_TREE}
  ${CMTK_INCLUDE_DIRS_SYSTEM}
  )

SUBDIRS(libs)

LINK_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR} ${DCMTK_LIBRARY_DIR}
  )

SET(CMTK_SYSTEM_LIBRARIES ${CMAKE_THREAD_LIBS})

#-----------------------------------------------------------------------------
# Dispatch the build into the proper subdirectories.

IF(BUILD_GUI)
  SUBDIRS(gui)
ENDIF(BUILD_GUI)    

SUBDIRS(scripts)

#-----------------------------------------------------------------------------
# Configure this tree

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.ctest.in ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.ctest @ONLY)

#-----------------------------------------------------------------------------
# Help other projects use CMTK.

# Copy the UseCMTK.cmake file to the binary tree for backward compatability.
CONFIGURE_FILE(${CMTK_SOURCE_DIR}/UseCMTK.cmake.in ${CMTK_BINARY_DIR}/UseCMTK.cmake COPYONLY IMMEDIATE)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/UseCMTK.cmake DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT libraries)

# Create the CMTKConfig.cmake file containing the CMTK configuration.
INCLUDE(${CMTK_SOURCE_DIR}/cmtkGenerateCMTKConfig.cmake)

#-----------------------------------------------------------------------------
# Add compiler flags CMTK needs to work on this platform.  This must be
# done AFTER the call to CMAKE_EXPORT_BUILD_SETTINGS.

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMTK_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMTK_REQUIRED_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMTK_REQUIRED_LINK_FLAGS}")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMTK_REQUIRED_LINK_FLAGS}")
SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${CMTK_REQUIRED_LINK_FLAGS}")

#-----------------------------------------------------------------------------
# Settings for rpath when building shared libraries

# Disable rpath altogether if we don't use shared libraries
IF(CMTK_BUILD_SHARED_LIBS)
  # use, i.e. don't skip the full RPATH for the build tree
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

  # the RPATH to be used when installing
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
ELSE(CMTK_BUILD_SHARED_LIBS)
  SET(CMAKE_SKIP_RPATH TRUE)
ENDIF(CMTK_BUILD_SHARED_LIBS)

# Save library dependencies.
EXPORT_LIBRARY_DEPENDENCIES(${CMTK_BINARY_DIR}/CMTKLibraryDepends.cmake)
INSTALL(FILES ${CMTK_BINARY_DIR}/CMTKLibraryDepends.cmake DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT libraries)

#-----------------------------------------------------------------------------
# After everything else is done, do we want to build man pages?

OPTION(BUILD_MANPAGES "Create application manpages" OFF)

