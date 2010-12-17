##
##  Copyright 2009-2010 SRI International
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

## Configure scripts for the build tree
SET(CMTK_BINARY_DIR_CONFIG ${CMAKE_BINARY_DIR}/bin)

SET(scripts utilities iterative_shape_averaging correct_nex_motion)
FOREACH(s ${scripts})
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${s}.sh.in ${CMAKE_BINARY_DIR}/bin/${s} @ONLY)
ENDFOREACH(s ${scripts})

## Configure scripts for the install tree
# Link directories.
SET(CMTK_BINARY_CONFIG "\${CMTK_INSTALL_PREFIX}${CMTK_INSTALL_BIN_DIR}")

FOREACH(s ${scripts})
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${s}.sh.in ${CMAKE_BINARY_DIR}/Utilities/${s} @ONLY)
  INSTALL(FILES ${CMAKE_BINARY_DIR}/Utilities/${s} DESTINATION bin COMPONENT runtime)
ENDFOREACH(s ${scripts})
