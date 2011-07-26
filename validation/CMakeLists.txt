##
##  Copyright 2009-2011 SRI International
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
SET(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW) ## libraries with full paths vs. search paths
ENDIF(COMMAND cmake_policy)
SET(CMAKE_OVERRIDE_COMPILER_MISMATCH 1)

PROJECT(CMTKvalidation)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMTKvalidation_SOURCE_DIR})

FIND_PACKAGE(CMTK)
IF(CMTK_FOUND)
  INCLUDE(${CMTK_USE_FILE})
ELSE(CMTK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build CMTKvalidation without CMTK.  Please set CMTK_DIR.")
ENDIF(CMTK_FOUND)

SUBDIRS(RIRE)
SUBDIRS(CURT)
