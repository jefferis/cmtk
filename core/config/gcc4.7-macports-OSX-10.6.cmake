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
##  $Revision: 3157 $
##
##  $LastChangedDate: 2011-04-18 13:37:45 -0700 (Mon, 18 Apr 2011) $
##
##  $LastChangedBy: torstenrohlfing $
##

SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
SET(CMAKE_INSTALL_PREFIX "/opt/local" CACHE PATH "Install prefix")

# General settings
SET(CMTK_SYSTEM_NAME "MacOSX-10.6-MacPorts" CACHE STRING "System name")
SET(CMTK_SYSTEM_PROCESSOR "x86_64" CACHE STRING "System processor")

# MacOS stuff
SET(CMAKE_OSX_SYSROOT "/" CACHE STRING "OS-X architecture")

# Select MacPorts compiler
SET(CMAKE_CXX_COMPILER "/opt/local/bin/g++-mp-4.7" CACHE FILEPATH "C++ compiler path")
SET(CMAKE_C_COMPILER "/opt/local/bin/gcc-mp-4.7" CACHE FILEPATH "C compiler path")

# Activate SSE support for floating point
SET(CMAKE_CXX_FLAGS "-march=core2 -mmmx -msse -msse2 -mfpmath=sse" CACHE STRING "C++ compiler flags")
SET(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS} CACHE STRING "C compiler flags")

# CMTK config settings
SET(CMTK_USE_OPENMP ON CACHE BOOL "Use OpenMP for parallelization" )
SET(CMTK_USE_GCD OFF CACHE BOOL "Use Grand Central Dispatch for SMP parallelism" )
SET(CMTK_USE_SQLITE OFF CACHE BOOL "Use SQLite database")
SET(CMTK_USE_QT ON CACHE BOOL "Use Qt toolkit")
SET(BUILD_GUI ON CACHE BOOL "Build GUI applications (requires Qt)")
SET(CMTK_USE_CUDA ON CACHE BOOL "Use CUDA for GPU acceleration" )
SET(CMTK_USE_LZMA OFF CACHE BOOL "Use LZMA library for decompression")

# Enable FFTW, since we require MacPorts for this build anyway
SET(CMTK_USE_FFTW ON CACHE BOOL "Use FFTW library (required for ADNI phantom detection)")

