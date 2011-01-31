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
##  $Revision: 2788 $
##
##  $LastChangedDate: 2011-01-20 14:02:10 -0800 (Thu, 20 Jan 2011) $
##
##  $LastChangedBy: torstenrohlfing $
##


FILE(APPEND "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "
CMAKE_OSX_ARCHITECTURES:STRING=i386
CMAKE_OSX_DEPLOYMENT_TARGET:STRING=10.4
CMAKE_OSX_SYSROOT:STRING=/Developer/SDKs/MacOSX10.4u.sdk
CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++-4.0
CMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc-4.0

CMAKE_CXX_FLAGS:STRING=-march=pentium4 -mmmx -msse -msse2 -mfpmath=sse -Wall -Wextra -Wno-deprecated -Wno-unknown-pragmas
CMAKE_C_FLAGS:STRING=-march=pentium4 -mmmx -msse -msse2 -mfpmath=sse

CMTK_USE_CUDA:BOOL=OFF
CMTK_USE_QT:BOOL=OFF
")
