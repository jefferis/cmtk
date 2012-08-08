##
##  Copyright 2010 Greg Jefferis
##
##  Copyright 2010-2012 SRI International
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
SET(CMAKE_INSTALL_PREFIX "/opt/local" CACHE PATH "Install prefix")

# General settings
SET(CMTK_SYSTEM_NAME "MacOSX-10.6" CACHE STRING "System name")
SET(CMTK_SYSTEM_PROCESSOR "x86_64" CACHE STRING "System processor")

# 64 bit for OS X >=10.6
SET(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "OS-X architecture")
SET(CMAKE_OSX_DEPLOYMENT_TARGET "10.6" CACHE STRING "OS-X target")
SET(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.6.sdk" CACHE PATH "OS-X SDK")

# Activate SSE support for floating point
SET(CMAKE_CXX_FLAGS "-march=core2 -mmmx -msse -msse2 -mfpmath=sse" CACHE STRING "C++ compiler flags")
SET(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS} CACHE STRING "C compiler flags")

# CMTK config settings
SET(CMTK_USE_SQLITE OFF CACHE BOOL "Use SQLite database")
SET(CMTK_USE_QT ON CACHE BOOL "Use Qt toolkit")
SET(BUILD_GUI ON CACHE BOOL "Build GUI applications (requires Qt)")
SET(CMTK_USE_CUDA ON CACHE BOOL "Use CUDA for GPU acceleration" )
SET(CMTK_USE_LZMA OFF CACHE BOOL "Use LZMA library for decompression")

# Disable OpenMP - broken on Mac
SET(CMTK_USE_OPENMP OFF CACHE BOOL "Use OpenMP for parallelization" )

# Disable FFTW, even if installed on our build system, since it requires MacPorts
SET(CMTK_USE_FFTW OFF CACHE BOOL "Use FFTW library (required for ADNI phantom detection)");
