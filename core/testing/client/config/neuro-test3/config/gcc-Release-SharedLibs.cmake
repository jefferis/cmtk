##
##  Copyright 1997-2011 Torsten Rohlfing
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

SET(CTEST_SITE "neuro-test3")
SET(CTEST_BUILD_NAME ${TEST_NAME})
SET(DART_TESTING_TIMEOUT 1800)

SET(CMTK_CTEST_ROOT "/Users/torsten/cmtk/testing/${TEST_NAME}")
SET(CTEST_SOURCE_DIRECTORY "${CMTK_CTEST_ROOT}/core")
SET(CTEST_BINARY_DIRECTORY "${CMTK_CTEST_ROOT}/build")
SET(CTEST_UPDATE_COMMAND "/usr/bin/svn")
SET(CTEST_BUILD_CONFIGURATION "Release")
SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")

IF(NOT EXISTS ${CTEST_SOURCE_DIRECTORY})
  SET(CTEST_CHECKOUT_COMMAND  "${CTEST_UPDATE_COMMAND} co https://www.nitrc.org:443/svn/cmtk/trunk/core/  ${CMTK_SOURCE_DIRECTORY}")
ENDIF(NOT EXISTS ${CTEST_SOURCE_DIRECTORY})
CTEST_EMPTY_BINARY_DIRECTORY(${CTEST_BINARY_DIRECTORY})

FILE(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "
BUILDNAME:STRING=Linux-gcc-x86_64-Release
BUILD_TESTING:BOOL=ON
BUILD_SHARED_LIBS:BOOL=ON
BUILD_GUI:BOOL=ON
BUILD_VALIDATION:BOOL=ON
CMAKE_BUILD_TYPE:STRING=Release
CMAKE_CXX_FLAGS_DEBUG:STRING=-g -DDEBUG

CMTK_USE_DCMTK:BOOL=ON
CMTK_USE_SQLITE:BOOL=ON
CMTK_BUILD_SQLITE:BOOL=ON
CMTK_USE_CUDA:BOOL=ON

CMTK_USE_QT:BOOL=ON
CMTK_BUILD_NRRD:BOOL=ON
")
