##
##  Copyright 1997-2009 Torsten Rohlfing
##
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

# turn off Windows security warnings for bundled code
IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

IF(CMTK_BUILD_ZLIB)
  SUBDIRS(zlib)
ENDIF(CMTK_BUILD_ZLIB)

IF(CMTK_BUILD_NRRD)
  SUBDIRS(NrrdIO)
ENDIF(CMTK_BUILD_NRRD)

IF(CMTK_BUILD_DCMTK)
  SUBDIRS(dcmtk)
ENDIF(CMTK_BUILD_DCMTK)

IF(CMTK_BUILD_MXML)
  SUBDIRS(mxml)
ENDIF(CMTK_BUILD_MXML)

IF(CMTK_BUILD_SQLITE)
  SUBDIRS(sqlite3)
ENDIF(CMTK_BUILD_SQLITE)
