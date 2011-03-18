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

SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")

# General settings
SET(CMTK_SYSTEM_NAME "MacOSX-10.4" CACHE STRING "System name")
SET(CMTK_SYSTEM_PROCESSOR "i686" CACHE STRING "System processor")

# 32 bit for OS X >=10.4
SET(CMAKE_OSX_ARCHITECTURES "i386" CACHE STRING "OS-X architectures")
SET(CMAKE_OSX_DEPLOYMENT_TARGET "10.4" CACHE STRING "OS-X target")
SET(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.5.sdk" CACHE PATH "OS-X SDK")

# Use GCC 4.0 for 10.4 SDK
SET(CMAKE_CXX_COMPILER "/usr/bin/g++-4.0" CACHE FILEPATH "C++ Compiler")
SET(CMAKE_C_COMPILER "/usr/bin/gcc-4.0" CACHE FILEPATH "C Compiler")

# Activate SSE support for floating point
SET(CMAKE_CXX_FLAGS "-m32 -march=pentium4 -mmmx -msse -msse2 -mfpmath=sse" CACHE STRING "C++ compiler flags")
SET(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS} CACHE STRING "C compiler flags")

# Disable Grand Central Dispatch as it is broken with C++ in 10.4 SDK
SET(CMTK_USE_GCD OFF CACHE BOOL "Use Grand Central Dispatch")
SET(CMTK_USE_OPENMP OFF CACHE BOOL "Use OpenMP for parallelization" )
