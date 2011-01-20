##
##  Copyright 2010 Greg Jefferis
##
##  Copyright 2010-2011 SRI International
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

# General settings
SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
SET(CMAKE_INSTALL_PREFIX "/opt/local" CACHE PATH "Install prefix")

# CMTK config settings
SET(CMTK_USE_SQLITE ON CACHE BOOL "Use SQLite database")
SET(CMTK_USE_QT OFF CACHE BOOL "Use Qt toolkit (does not work for 10.4 right now)")
SET(BUILD_GUI OFF CACHE BOOL "Build GUI applications (requires Qt)")
SET(CMTK_USE_LZMA OFF CACHE BOOL "Use LZMA library for decompression")

# 32 bit for OS X >=10.4
SET(CMAKE_OSX_ARCHITECTURES "i386" CACHE STRING "OS-X architectures")
SET(CMAKE_OSX_DEPLOYMENT_TARGET "10.4" CACHE STRING "OS-X target")
SET(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.4u.sdk" CACHE PATH "OS-X SDK")

# Activate SSE support for floating point
SET(CMAKE_CXX_FLAGS "-march=pentium4 -mmmx -msse -msse2 -mfpmath=sse" CACHE STRING "C compiler flags")
SET(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS} CACHE STRING "C++ compiler flags")
