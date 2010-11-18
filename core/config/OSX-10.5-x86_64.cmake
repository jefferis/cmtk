##
##  Copyright 2010 Greg Jefferis
##
##  Copyright 2010 SRI International
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

# 64 bit for OS X >=10.5
SET(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "OS-X architecture")
SET(CMAKE_OSX_DEPLOYMENT_TARGET "10.5" CACHE STRING "OS-X target")
SET(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.5.sdk" CACHE PATH "OS-X SDK")

# Activate SSE support for floating point
SET(CMAKE_CXX_FLAGS "-march=nocona -mmmx -msse -msse2 -mfpmath=sse" CACHE STRING "C compiler flags")
SET(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS} CACHE STRING "C++ compiler flags")
