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

SET(TEST_NAME i686-Release)
SET(CTEST_SITE "neuro-test0")
SET(CTEST_BUILD_NAME ${TEST_NAME})
SET(DART_TESTING_TIMEOUT 1800)

SET(CMTK_CTEST_ROOT "/home/testrunner/nitrc/${TEST_NAME}")
SET(CTEST_PROJECT_NAME "CMTK")
SET(CTEST_SOURCE_DIRECTORY "${CMTK_CTEST_ROOT}/core")
SET(CTEST_BINARY_DIRECTORY "${CMTK_CTEST_ROOT}/build")
SET(CTEST_UPDATE_COMMAND "/usr/bin/svn")
SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")
SET(CTEST_BUILD_CONFIGURATION "Release")

IF(NOT EXISTS ${CTEST_SOURCE_DIRECTORY})
  SET(CTEST_CHECKOUT_COMMAND  "${CTEST_UPDATE_COMMAND} co https://www.nitrc.org:443/svn/cmtk/trunk/core/  ${CMTK_SOURCE_DIRECTORY}")
ENDIF(NOT EXISTS ${CTEST_SOURCE_DIRECTORY})
CTEST_EMPTY_BINARY_DIRECTORY(${CTEST_BINARY_DIRECTORY})

FILE(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "
BUILDNAME:STRING=Linux-c++-i686-Release
BUILD_TESTING:BOOL=ON
CMTK_DATA_ROOT:PATH=/home/testrunner/nitrc/data
BUILD_FUSION:BOOL=ON
BUILD_GUI:BOOL=ON
BUILD_VALIDATION:BOOL=ON
CMAKE_BUILD_TYPE:STRING=Release
CMAKE_CXX_FLAGS:STRING=-m32 -march=pentium4 -mmmx -msse -msse2 -mfpmath=sse -Wno-deprecated -Wno-unknown-pragmas
CMAKE_CXX_FLAGS_DEBUG:STRING=-g -DDEBUG
CMAKE_C_FLAGS:STRING=-m32 -march=pentium4 -mmmx -msse -msse2 -mfpmath=sse
DCMTK_INCLUDE_DIR:PATH=/usr/include
DCMTK_LIBRARY_DIR:PATH=/usr/lib
CMTK_USE_DCMTK:BOOL=OFF
CMTK_USE_QT:BOOL=ON
QT_INCLUDE_DIR:PATH=/usr/include
QT_LIBRARY_DIR:PATH=/usr/lib/
CMTK_BUILD_NRRD:BOOL=ON
")

CTEST_START(Continuous)
CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
IF(${res} GREATER 0)
	CTEST_CONFIGURE(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
	CTEST_BUILD(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
	CTEST_TEST(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
	CTEST_MEMCHECK(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
	CTEST_COVERAGE(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
	CTEST_SUBMIT(RETURN_VALUE res)
ENDIF(${res} GREATER 0)
