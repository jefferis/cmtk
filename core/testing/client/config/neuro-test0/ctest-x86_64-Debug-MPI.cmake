##
##  Copyright 1997-2009 Torsten Rohlfing
##  Copyright 2004-2009 SRI International
##
##  This file is part of the Computational Morphometry Toolkit.
##
##  http:##www.nitrc.org/projects/cmtk/
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
##  <http:##www.gnu.org/licenses/>.
##
##  $Revision$
##
##  $LastChangedDate$
##
##  $LastChangedBy$
##

SET (TEST_NAME x86_64-Debug-MPI)
SET(CTEST_SITE "neuro-test0")
SET(CTEST_BUILD_NAME ${TEST_NAME})
SET(DART_TESTING_TIMEOUT 1800)

SET (CTEST_SOURCE_DIRECTORY "/home/testrunner/source/${TEST_NAME}")
SET (CTEST_BINARY_DIRECTORY "/home/testrunner/builds/${TEST_NAME}")
SET(CTEST_UPDATE_COMMAND                "/usr/bin/svn")
SET(CTEST_BUILD_CONFIGURATION           "Debug")

IF(NOT EXISTS ${CTEST_SOURCE_DIRECTORY})
  SET(CTEST_CHECKOUT_COMMAND  "${CTEST_UPDATE_COMMAND} co https://neuro.sri.com/repos/igs/trunk ${CTEST_SOURCE_DIRECTORY}")
ENDIF(NOT EXISTS ${CTEST_SOURCE_DIRECTORY})

SET(CTEST_COVERAGE_COMMAND              "/usr/bin/gcov")
SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")

CTEST_EMPTY_BINARY_DIRECTORY(${CTEST_BINARY_DIRECTORY})

FILE(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "
BUILDNAME:STRING=Linux-mpic++-x86_64-Debug

CMAKE_C_COMPILER:PATH=/usr/bin/mpicc
CMAKE_CXX_COMPILER:PATH=/usr/bin/mpic++

BUILD_FUSION:BOOL=OFF
BUILD_GUI:BOOL=OFF
CMAKE_BUILD_TYPE:STRING=Debug
CMAKE_CXX_FLAGS:STRING=-m64 -march=nocona -mmmx -msse -msse2 -mfpmath=sse -Wno-deprecated -Wno-unknown-pragmas
CMAKE_CXX_FLAGS_DEBUG:STRING=-g -DDEBUG -fprofile-arcs -ftest-coverage
CMAKE_C_FLAGS:STRING=-m64 -march=nocona -mmmx -msse -msse2 -mfpmath=sse
CMAKE_C_FLAGS_DEBUG:STRING=-g -DDEBUG -fprofile-arcs -ftest-coverage
DCMTK_INCLUDE_DIR:PATH=/usr/include
DCMTK_LIBRARY_DIR:PATH=/usr/lib64
CMTK_USE_MPI:BOOL=ON
CMTK_USE_DCMTK:BOOL=ON
CMTK_USE_QT:BOOL=OFF
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
