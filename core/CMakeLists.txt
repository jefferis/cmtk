##
##  Copyright 1997-2009 Torsten Rohlfing
##  Copyright 2004-2009 SRI International
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW) ## libraries with full paths vs. search paths
  cmake_policy(SET CMP0005 OLD) ## escaping of strings in ADD_DEFINITIONS
ENDIF(COMMAND cmake_policy)
SET(CMAKE_OVERRIDE_COMPILER_MISMATCH 1)

PROJECT(CMTK)
SET(CMTK_VERSION_MAJOR "1")
SET(CMTK_VERSION_MINOR "1")
SET(CMTK_VERSION_PATCH "1-svn")

SET(CMTK_VERSION_STRING "${CMTK_VERSION_MAJOR}.${CMTK_VERSION_MINOR}.${CMTK_VERSION_PATCH}")
ADD_DEFINITIONS(-DCMTK_VERSION=\\\"${CMTK_VERSION_STRING}\\\")

#-----------------------------------------------------------------------------
# Test for some required system information.
INCLUDE (${CMAKE_ROOT}/Modules/CMakeDetermineSystem.cmake)
INCLUDE (${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeSystem.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityC.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityCXX.cmake)

ADD_DEFINITIONS(-DCMTK_INSTANTIATE_TEMPLATES)
#-----------------------------------------------------------------------------
# Check for system headers.
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES (dirent.h HAVE_DIRENT_H)
CHECK_INCLUDE_FILES (execinfo.h HAVE_EXECINFO_H)
CHECK_INCLUDE_FILES (fcntl.h HAVE_FCNTL_H)
CHECK_INCLUDE_FILES (float.h HAVE_FLOAT_H)
CHECK_INCLUDE_FILES (ieeefp.h HAVE_IEEEFP_H)
CHECK_INCLUDE_FILES (limits.h HAVE_LIMITS_H)
CHECK_INCLUDE_FILES (malloc.h HAVE_MALLOC_H)
CHECK_INCLUDE_FILES (pthread.h HAVE_PTHREAD_H)
CHECK_INCLUDE_FILES (signal.h HAVE_SIGNAL_H)
CHECK_INCLUDE_FILES (stdarg.h HAVE_STDARG_H)
CHECK_INCLUDE_FILES (stddef.h HAVE_STDDEF_H)
CHECK_INCLUDE_FILES (stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILES (stdio.h HAVE_STDIO_H)
CHECK_INCLUDE_FILES (stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILES (string.h HAVE_STRING_H)
CHECK_INCLUDE_FILES (time.h HAVE_TIME_H)
CHECK_INCLUDE_FILES (ucontext.h HAVE_UCONTEXT_H)
CHECK_INCLUDE_FILES (unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILES (values.h HAVE_VALUES_H)
CHECK_INCLUDE_FILES (varargs.h HAVE_VARARGS_H)

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
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles/TestQnanhibit.cmake)
TEST_QNANHIBIT(QNANHIBIT_VALUE ${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles)
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

OPTION(CMTK_BUILD_UNSTABLE "Build 'unstable' library and applications" OFF)
MARK_AS_ADVANCED(CMTK_BUILD_UNSTABLE)
IF(CMTK_BUILD_UNSTABLE)
  ADD_DEFINITIONS(-DCMTK_BUILD_UNSTABLE)
ENDIF(CMTK_BUILD_UNSTABLE)

OPTION(CMTK_BUILD_STACKTRACE "Build applications with stack trace printing in case of crash (increases build time)" ON)
MARK_AS_ADVANCED(CMTK_BUILD_STACKTRACE)

OPTION(CMTK_SINGLE_COMMAND_BINARY "Build using a single wrapper binary for all command line tools to avoid name collisions with other packages" OFF)
MARK_AS_ADVANCED(CMTK_SINGLE_COMMAND_BINARY)
IF(CMTK_SINGLE_COMMAND_BINARY)
  ADD_DEFINITIONS(-DCMTK_SINGLE_COMMAND_BINARY)
ENDIF(CMTK_SINGLE_COMMAND_BINARY)

#-----------------------------------------------------------------------------
# Configure Dart testing support.  This should be done before any
# MESSAGE(FATAL_ERROR ...) commands are invoked.

INCLUDE(${CMAKE_ROOT}/Modules/Dart.cmake)
MARK_AS_ADVANCED(TCL_TCLSH DART_ROOT)
ENABLE_TESTING()
IF(BUILD_TESTING)
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)

  OPTION(CMTK_TESTING_MEMORYCHECK "Test using memory checking" OFF)
  MARK_AS_ADVANCED(CMTK_TESTING_MEMORYCHECK)

  IF(NOT CMTK_DATA_ROOT)
    SET(CMTK_DATA_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../data CACHE PATH "Root directory of the CMTK data tree")
  ENDIF(NOT CMTK_DATA_ROOT)
  MARK_AS_ADVANCED(CMTK_DATA_ROOT)
ENDIF(BUILD_TESTING)

SET(CMTK_ROOT_PATH_SRI24 "CMTK_ROOT_PATH_SRI24-NOTFOUND" CACHE PATH "Root path of the SRI24 atlas data (usually ends with \"sri24\")" )

#-----------------------------------------------------------------------------
# Configure Doxygen support

FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
  OPTION(BUILD_DOCUMENTATION "Build project documentation using Doxygen" OFF)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  IF(BUILD_DOCUMENTATION)
    ADD_CUSTOM_TARGET(doc COMMAND doxygen DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    INSTALL(DIRECTORY html/ DESTINATION share/doc/cmtk/)
  ENDIF(BUILD_DOCUMENTATION)
ENDIF(DOXYGEN_FOUND)

#-----------------------------------------------------------------------------
# Add an option to use double-precision floats for coordinates.

OPTION(CMTK_COORDINATES_DOUBLE "Use double precision floats for coordinates (turn off to save memory)" ON)
IF(CMTK_COORDINATES_DOUBLE)
  ADD_DEFINITIONS(-DCMTK_COORDINATES_DOUBLE)
ELSE(CMTK_COORDINATES_DOUBLE)
  ADD_DEFINITIONS(-DCMTK_COORDINATES_FLOAT)
ENDIF(CMTK_COORDINATES_DOUBLE)
MARK_AS_ADVANCED(CMTK_COORDINATES_DOUBLE) 

OPTION(CMTK_DATA_DOUBLE "Use double precision floats for data exchange between typed arrays" ON)
IF(CMTK_DATA_DOUBLE)
  ADD_DEFINITIONS(-DCMTK_DATA_DOUBLE)
ELSE(CMTK_DATA_DOUBLE)
  ADD_DEFINITIONS(-DCMTK_DATA_FLOAT)
ENDIF(CMTK_DATA_DOUBLE)
MARK_AS_ADVANCED(CMTK_DATA_DOUBLE) 

OPTION(CMTK_NUMERICS_DOUBLE "Use double precision floats for numerical algorithms" ON)
IF(CMTK_NUMERICS_DOUBLE)
  ADD_DEFINITIONS(-DCMTK_NUMERICS_DOUBLE)
ELSE(CMTK_NUMERICS_DOUBLE)
  ADD_DEFINITIONS(-DCMTK_NUMERICS_FLOAT)
ENDIF(CMTK_NUMERICS_DOUBLE)
MARK_AS_ADVANCED(CMTK_NUMERICS_DOUBLE)

#-----------------------------------------------------------------------------
# Add an option to use or not use DCMTK

OPTION(CMTK_USE_DCMTK "Use DCMTK (for DICOM import)" OFF)

#-----------------------------------------------------------------------------
# Add an option to use or not use DCMTK

OPTION(CMTK_USE_BUNDLED_DCMTK "Accept the bundled version of DCMTK" OFF)

IF(CMTK_USE_DCMTK)
  IF(CMTK_USE_BUNDLED_DCMTK)
    MESSAGE( STATUS "Building bundled DCMTK library" )
    SET(CMTK_BUILD_DCMTK 1)
    SET(HAVE_DCMTK 1)
    SET(DCMTK_SOURCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/Utilities/dcmtk)
    SET(DCMTK_BINARY_ROOT ${CMAKE_CURRENT_BINARY_DIR}/Utilities/dcmtk)
    SET(DCMTK_INCLUDE_DIR "${DCMTK_SOURCE_ROOT}")
    SET(DCMTK_INCLUDE_DIRS "${DCMTK_BINARY_ROOT};${DCMTK_INCLUDE_DIR};${DCMTK_INCLUDE_DIR}/dcmtk/dcmjpeg;${DCMTK_INCLUDE_DIR}/dcmtk/dcmimage;${DCMTK_INCLUDE_DIR}/dcmtk/dcmimgle;${DCMTK_INCLUDE_DIR}/dcmtk/dcmdata;${DCMTK_INCLUDE_DIR}/dcmtk/ofstd")
    SET(DCMTK_LIBRARIES "D4CMTK")
    SET(DCMTK_DICTIONARY_DIR ${DCMTK_SOURCE_ROOT})

    ADD_DEFINITIONS(-DCMTK_HAVE_DCMTK -DCMAKE_DCMDICTPATH=\\\"${DCMTK_DICTIONARY_DIR}\\\")

    IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
      SET(DCMTK_LIBRARIES_EXTRA "wsock32;netapi32" CACHE STRING "Extra libraries needed for DCMTK")
    ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
      ADD_DEFINITIONS(-DHAVE_CONFIG_H)
    ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

    SET(DCM_DICT_DEFAULT_PATH "${DCMTK_DICTIONARY_DIR}")
  ELSE (CMTK_USE_BUNDLED_DCMTK)
    IF(NOT DCMTK_DIR)
      SET(DCMTK_DIR "NOTFOUND")
    ENDIF(NOT DCMTK_DIR)
    
    FIND_PATH(DCMTK_LIBRARY_DIR dicom.dic ${DCMTK_DIR}/lib /usr/local/lib/ /usr/lib/)
    IF(DCMTK_LIBRARY_DIR)
      SET(DCMTK_DICTIONARY_DIR ${DCMTK_LIBRARY_DIR})

      # Install dicom.dic from dcmtk so it's in packages
      IF(CMTK_USE_DCMTK)
	INSTALL(FILES ${DCMTK_DICTIONARY_DIR}/dicom.dic DESTINATION lib)
      ENDIF(CMTK_USE_DCMTK)

      SET(DCMTK_LIBRARIES "dcmjpeg;ijg8;ijg12;ijg16;dcmimage;dcmimgle;dcmdata;ofstd")
      MARK_AS_ADVANCED(DCMTK_LIBRARY_DIR)
    ELSE(DCMTK_LIBRARY_DIR)
      MESSAGE(FATAL_ERROR "Please set DCMTK_LIBRARY_DIR")		
    ENDIF(DCMTK_LIBRARY_DIR)
    
    FIND_PATH(DCMTK_INCLUDE_DIR dcmtk ${DCMTK_DIR}/include /usr/local/include /usr/include )
    IF(DCMTK_INCLUDE_DIR)
      SET(DCMTK_INCLUDE_DIRS "${DCMTK_INCLUDE_DIR};${DCMTK_INCLUDE_DIR}/dcmtk/dcmjpeg;${DCMTK_INCLUDE_DIR}/dcmtk/ijg8;${DCMTK_INCLUDE_DIR}/dcmtk/ijg12;${DCMTK_INCLUDE_DIR}/dcmtk/ijg16;${DCMTK_INCLUDE_DIR}/dcmtk/dcmimage;${DCMTK_INCLUDE_DIR}/dcmtk/dcmimgle;${DCMTK_INCLUDE_DIR}/dcmtk/dcmdata;${DCMTK_INCLUDE_DIR}/dcmtk/ofstd")
      MARK_AS_ADVANCED(DCMTK_INCLUDE_DIR)
    ELSE(DCMTK_INCLUDE_DIR)
      MESSAGE(FATAL_ERROR "Please set DCMTK_INCLUDE_DIR")		
    ENDIF(DCMTK_INCLUDE_DIR)
    
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
      SET(DCMTK_LIBRARIES_EXTRA "wsock32;netapi32" CACHE STRING "Extra libraries needed for DCMTK")
    ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
      SET(DCMTK_LIBRARIES_EXTRA "" CACHE STRING "Extra libraries needed for DCMTK")
    ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    
    MARK_AS_ADVANCED(DCMTK_LIBRARIES_EXTRA)
    
    IF(DCMTK_LIBRARY_DIR AND DCMTK_INCLUDE_DIR)
      ADD_DEFINITIONS(-DCMTK_HAVE_DCMTK -DCMAKE_DCMDICTPATH=\\\"${DCMTK_LIBRARY_DIR}/dicom.dic\\\")
      IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        ADD_DEFINITIONS(-DHAVE_CONFIG_H)
      ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    ENDIF(DCMTK_LIBRARY_DIR AND DCMTK_INCLUDE_DIR)
  ENDIF (CMTK_USE_BUNDLED_DCMTK)
ENDIF(CMTK_USE_DCMTK)

#-----------------------------------------------------------------------------
# Add an option to use bundled zlib.

# Check for useful libraries
Find_Package(ZLIB)
IF(NOT ZLIB_FOUND)
  MESSAGE( STATUS "Building bundled zlib library" )
  SET(CMTK_BUILD_ZLIB 1)
  SET(ZLIB_FOUND 1)
  SET(ZLIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utilities/zlib CACHE PATH "Directory containing zlib.h header file" FORCE)
  SET(ZLIB_LIBRARIES "cmtkZlib" CACHE STRING "Libraries to link for zlib support")
ENDIF(NOT ZLIB_FOUND)

#-----------------------------------------------------------------------------
# Add an option to use or not use TEEM (for NrrdIO support)

OPTION(CMTK_BUILD_NRRD "Build NrrdIO library for NRRD import/export" ON)

IF(CMTK_BUILD_NRRD)
  ADD_DEFINITIONS(-DCMTK_BUILD_NRRD -DTEEM_ZLIB)
  SET(NRRD_LIBRARIES "NrrdIO")
  SET(NRRD_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utilities/NrrdIO)
ENDIF(CMTK_BUILD_NRRD)

#-----------------------------------------------------------------------------
# Configure bundled utilities.

ADD_SUBDIRECTORY(Utilities)

SET(MXML_LIBRARIES "mxml")
SET(MXML_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Utilities/mxml")


#-----------------------------------------------------------------------------
# Add an option to use or not use VTK

OPTION(CMTK_USE_VTK "Use VTK (The Visualization Toolkit) (some optional GUI applications need this)" OFF)
IF(CMTK_USE_VTK)
  FIND_PACKAGE(VTK)
  IF (VTK_FOUND)
    INCLUDE (${VTK_USE_FILE})

    ADD_DEFINITIONS(-DCMTK_HAVE_VTK)
    SET(CMTK_VTK_LIBS vtkHybrid vtkGraphics vtkCommon)
    # Test for VTK using std libs
    IF(WIN32) # It seems that the following message is only relavant for MSVC 6 installs
      IF(NOT VTK_USE_ANSI_STDLIB)
        MESSAGE("Warning.  Your VTK was not built with the VTK_USE_ANSI_STDLIB "
          "option ON.  Link errors may occur.  Consider re-building VTK "
          "with this option ON.  For MSVC 6, you MUST turn on the VTK "
          "option. CMTK will not link otherwise.")
      ENDIF(NOT VTK_USE_ANSI_STDLIB)
    ENDIF(WIN32)
  ENDIF (VTK_FOUND)
ENDIF(CMTK_USE_VTK)

#-----------------------------------------------------------------------------
# Add an option to enable or disable QT support and GUI applications

IF(NOT CMTK_USE_MPI)
  OPTION(CMTK_USE_QT "Use Qt (GUI applications need this)" OFF)
  IF(CMTK_USE_QT)
    SET (QT_MT_REQUIRED true)
    SET (QT_MIN_VERSION "4.3.0")
    FIND_PACKAGE(Qt4 COMPONENTS QtCore QtGui Qt3Support REQUIRED)
    IF (QT_FOUND)
      MESSAGE("Qt4 Found OK (${qt_version_str})")

      SET(QT_USE_QT3SUPPORT 1)
      INCLUDE(${QT_USE_FILE})
    ELSE(QT_FOUND)
      MESSAGE(FATAL_ERROR "No Qt4")
    ENDIF(QT_FOUND)

    ADD_DEFINITIONS(-DQT3_SUPPORT)

  ENDIF(CMTK_USE_QT)
ENDIF(NOT CMTK_USE_MPI)

#-----------------------------------------------------------------------------
# Add an option to use or not use AMD Core Math Library

OPTION(CMTK_USE_ACML "Use AMD Core Math Library (ACML) for better math performance" OFF)
IF(CMTK_USE_ACML)
  SET(ACML_FOUND 1)
  FIND_PATH(ACML_INCLUDE_DIR acml.h ${INCLUDE_DIRS})
  IF(NOT EXISTS ${ACML_INCLUDE_DIR}/acml_mv.h)
    SET(ACML_INCLUDE_DIR "NOTFOUND" CACHE PATH "Include directory of ACML library" )
    SET(ACML_FOUND 0)
  ENDIF(NOT EXISTS ${ACML_INCLUDE_DIR}/acml_mv.h)

  FIND_PATH(ACML_LIBRARY_DIR libacml.a ${LIBRARY_DIRS})
  IF(NOT EXISTS ${ACML_LIBRARY_DIR}/libacml_mv.a)
    SET(ACML_LIBRARY_DIR "NOTFOUND" CACHE PATH "Include directory of ACML library" )
    SET(ACML_FOUND 0)
  ENDIF(NOT EXISTS ${ACML_LIBRARY_DIR}/libacml_mv.a)

  IF(ACML_FOUND)
    ADD_DEFINITIONS(-DCMTK_HAVE_ACML)
    SET(ACML_LIBRARIES "acml_mv")
    MARK_AS_ADVANCED(ACML_INCLUDE_DIR ACML_LIBRARY_DIR)
  ENDIF(ACML_FOUND)

ENDIF(CMTK_USE_ACML)

#-----------------------------------------------------------------------------
# Add an option to use or not to use MPI distributed-memory parallelism

OPTION(CMTK_USE_MPI "Use MPI for distributed-memory parallelism" OFF)
IF(CMTK_USE_MPI)
  ADD_DEFINITIONS(-DCMTK_BUILD_MPI)
ENDIF(CMTK_USE_MPI)

OPTION(CMTK_BUILD_SMP "Enable shared-memory parallelism (requires at least POSIX Threads, OpenMP, or Windows threads)" ON)
IF(CMTK_BUILD_SMP)
  ADD_DEFINITIONS(-DCMTK_BUILD_SMP)
ENDIF(CMTK_BUILD_SMP)


#-----------------------------------------------------------------------------
# Add an option to enable or disable POSIX threads if pthread.h header exists

IF(HAVE_PTHREAD_H)
  OPTION(CMTK_USE_PTHREADS "Use POSIX Threads for SMP parallelism" ON)
  MARK_AS_ADVANCED(CMTK_USE_PTHREADS)

  IF(NOT CMTK_BUILD_SMP)
    SET(CMTK_USE_PTHREADS OFF)
  ENDIF(NOT CMTK_BUILD_SMP)
ENDIF(HAVE_PTHREAD_H)

#-----------------------------------------------------------------------------
# for the gnu compiler a -D_PTHREADS is needed on sun
# for the native compiler a -mt flag is needed on the sun
IF(CMTK_USE_PTHREADS)

  ADD_DEFINITIONS(-DCMTK_USE_THREADS)

  IF(CMAKE_SYSTEM MATCHES "SunOS.*")
    IF(CMAKE_COMPILER_IS_GNUCXX)
      SET(CMTK_REQUIRED_CXX_FLAGS "${CMTK_REQUIRED_CXX_FLAGS} -D_PTHREADS")
    ELSE(CMAKE_COMPILER_IS_GNUCXX)
      SET(CMTK_REQUIRED_CXX_FLAGS "${CMTK_REQUIRED_CXX_FLAGS} -mt")
      SET(CMTK_REQUIRED_C_FLAGS "${CMTK_REQUIRED_C_FLAGS} -mt")
    ENDIF(CMAKE_COMPILER_IS_GNUCXX)
  ENDIF(CMAKE_SYSTEM MATCHES "SunOS.*")
ENDIF(CMTK_USE_PTHREADS)

# mingw thread support
IF(MINGW)
  SET(CMTK_REQUIRED_CXX_FLAGS "${CMTK_REQUIRED_CXX_FLAGS} -mthreads")
  SET(CMTK_REQUIRED_C_FLAGS "${CMTK_REQUIRED_C_FLAGS} -mthreads")
  SET(CMTK_REQUIRED_LINK_FLAGS "${CMTK_REQUIRED_LINK_FLAGS} -mthreads")
ENDIF(MINGW)

#-----------------------------------------------------------------------------
# Check for, and configure, OpenMP

IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION} STRLESS 2.6.2)
  MESSAGE("Warning.  We need at least CMake 2.6.2 to check for OpenMP support. You have version ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}")
ELSE (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION} STRLESS 2.6.2)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    OPTION(CMTK_USE_OPENMP "Use OpenMP for SMP parallelism" ON)
    MARK_AS_ADVANCED(CMTK_USE_OPENMP)

    IF(NOT CMTK_BUILD_SMP)
      SET(CMTK_USE_OPENMP OFF)
    ENDIF(NOT CMTK_BUILD_SMP)

    IF(CMTK_USE_OPENMP)
      SET(CMTK_REQUIRED_CXX_FLAGS "${CMTK_REQUIRED_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      SET(CMTK_REQUIRED_LINK_FLAGS "${CMTK_REQUIRED_LINK_FLAGS} ${OpenMP_C_FLAGS}")
      ADD_DEFINITIONS(-DCMTK_USE_OPENMP)
    ENDIF(CMTK_USE_OPENMP)
  ENDIF(OPENMP_FOUND)
ENDIF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION} STRLESS 2.6.2)

#-----------------------------------------------------------------------------
# Check whether compiler supports variable-size automatic arrays, i.e.,
# things like "const int size = 10; const char string[size];" rather than
# "const char string[10];"

INCLUDE(${CMAKE_ROOT}/Modules/CheckCXXSourceCompiles.cmake)
CHECK_CXX_SOURCE_COMPILES("int main() { int size = 10; char string[size]; }" COMPILER_VAR_AUTO_ARRAYSIZE)
IF(COMPILER_VAR_AUTO_ARRAYSIZE)
  ADD_DEFINITIONS(-DCMTK_VAR_AUTO_ARRAYSIZE)
ENDIF(COMPILER_VAR_AUTO_ARRAYSIZE)

CHECK_CXX_SOURCE_COMPILES("template<int N> class C { public: static const float nn = 1.0 / N; }; int main() { C<10> c; return C<10>::nn; }; " COMPILER_INCLASS_MEMBER_INITIALIZATION)
IF(NOT COMPILER_INCLASS_MEMBER_INITIALIZATION)
  ADD_DEFINITIONS(-DCMTK_NO_INCLASS_MEMBER_INITIALIZATION)
ENDIF(NOT COMPILER_INCLASS_MEMBER_INITIALIZATION)

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
  SET(CMTK_INSTALL_BIN_DIR /bin CACHE PATH "Binary files will be installed here")
ENDIF(NOT CMTK_INSTALL_BIN_DIR)

#-----------------------------------------------------------------------------
# Set up CPack support, to be able to distribute binary packages.
#
INCLUDE(InstallRequiredSystemLibraries)
IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  SET(CPACK_GENERATOR "ZIP;NSIS")
ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CPACK_GENERATOR "TGZ;RPM")
  ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CPACK_GENERATOR "TGZ")
  ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CMTK -- The Computational Morphometry Toolkit")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_PACKAGE_VENDOR "SRI International - Neuroscience Program")
SET(CPACK_PACKAGE_VERSION_MAJOR "${CMTK_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${CMTK_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${CMTK_VERSION_PATCH}")
SET(CPACK_PACKAGE_FILE_NAME "CMTK-${CMTK_VERSION_MAJOR}.${CMTK_VERSION_MINOR}.${CMTK_VERSION_PATCH}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMTK-${CMTK_VERSION_MAJOR}.${CMTK_VERSION_MINOR}.${CMTK_VERSION_PATCH}-${CMAKE_SYSTEM_NAME}/${CMAKE_SYSTEM_PROCESSOR}")
INCLUDE(CPack)

INSTALL(FILES ${CPACK_RESOURCE_FILE_LICENSE} ${CPACK_PACKAGE_DESCRIPTION_FILE} DESTINATION share/doc/CMTK)

#-----------------------------------------------------------------------------
# Dispatch the build into the proper subdirectories.

#-----------------------------------------------------------------------------
# Create the list of include directories needed for header files.
INCLUDE(${CMTK_SOURCE_DIR}/cmtkIncludeDirectories.cmake)

# This should be the only INCLUDE_DIRECTORIES command in the entire
# tree.
INCLUDE_DIRECTORIES(
  ${CMTK_INCLUDE_DIRS_BUILD_TREE}
  ${CMTK_INCLUDE_DIRS_BUILD_TREE_CXX}
  ${ACML_INCLUDE_DIR} ${DCMTK_INCLUDE_DIRS}
  ${QT_INCLUDE_DIR} ${QT_INCLUDE_PATH}
  ${NRRD_INCLUDE_DIR}
  ${MXML_INCLUDE_DIR}
  ${ZLIB_INCLUDE_DIR}
  ${CMTK_INCLUDE_DIRS_SYSTEM}
  )

SUBDIRS(libs)

OPTION(BUILD_APPS "Build command line applications" ON)
IF(BUILD_APPS)
  SUBDIRS(apps)
ENDIF(BUILD_APPS)

OPTION(BUILD_VALIDATION "Build tools for algorithm validation" OFF)
IF(BUILD_VALIDATION)
  SUBDIRS(validation)
ENDIF(BUILD_VALIDATION)

IF (BUILD_GUI)
  SUBDIRS(gui)
ENDIF (BUILD_GUI)    

LINK_DIRECTORIES(
  ${ACML_LIBRARY_DIR} ${DCMTK_LIBRARY_DIR}
  )

LINK_LIBRARIES(
  ${DCMTK_LIBRARIES} ${DCMTK_LIBRARIES_EXTRA} ${NRRD_LIBRARIES} ${MXML_LIBRARIES} ${ACML_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS}
  )

#-----------------------------------------------------------------------------
# Add options to enable or disable GUI applications when building with Qt
# support

IF(NOT CMTK_USE_MPI)
  IF(CMTK_USE_QT)
    OPTION(BUILD_GUI "Build GUI applications" ON)
    OPTION(BUILD_FUSION "Build image fusion application" ON)
  ENDIF(CMTK_USE_QT)
ENDIF(NOT CMTK_USE_MPI)

OPTION(BUILD_TESTING "Build test code" ON)
IF (BUILD_TESTING)
  SUBDIRS(testing)
ENDIF (BUILD_TESTING)

#-----------------------------------------------------------------------------
# Configure this tree

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/cmtkconfig.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/cmtkconfig.h)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.ctest.in ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.ctest @ONLY)

#-----------------------------------------------------------------------------
# Help other projects use CMTK.

# Copy the UseCMTK.cmake file to the binary tree for backward compatability.
CONFIGURE_FILE(${CMTK_SOURCE_DIR}/UseCMTK.cmake.in ${CMTK_BINARY_DIR}/UseCMTK.cmake COPYONLY IMMEDIATE)

# Save the compiler settings so another project can import them.
INCLUDE(${CMAKE_ROOT}/Modules/CMakeExportBuildSettings.cmake)
CMAKE_EXPORT_BUILD_SETTINGS(${CMTK_BINARY_DIR}/CMTKBuildSettings.cmake)

# Save library dependencies.
EXPORT_LIBRARY_DEPENDENCIES(${CMTK_BINARY_DIR}/CMTKLibraryDepends.cmake)

# Create the CMTKConfig.cmake file containing the CMTK configuration.
INCLUDE (${CMTK_SOURCE_DIR}/cmtkGenerateCMTKConfig.cmake)

#-----------------------------------------------------------------------------
# Add compiler flags CMTK needs to work on this platform.  This must be
# done AFTER the call to CMAKE_EXPORT_BUILD_SETTINGS.

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMTK_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMTK_REQUIRED_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMTK_REQUIRED_LINK_FLAGS}")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMTK_REQUIRED_LINK_FLAGS}")
SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${CMTK_REQUIRED_LINK_FLAGS}")

# Disable rpath altogether if we don't use shared libraries
IF(${CMTK_BUILD_SHARED_LIBS})
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
ELSE(${CMTK_BUILD_SHARED_LIBS})
  SET(CMAKE_SKIP_RPATH TRUE)
ENDIF(${CMTK_BUILD_SHARED_LIBS})
