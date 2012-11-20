##
##  Copyright 2012 SRI International
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

SET(TEMPLATES structural_pipeline_index.html structural.tex)

FOREACH(t ${TEMPLATES})
  ## Configure for the build tree
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${t}.in ${CMAKE_CURRENT_BINARY_DIR}/${t} @ONLY)
  
  ## Configure for the install tree
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${t}.in ${CMAKE_BINARY_DIR}/Install/ncanda/templates/${t} @ONLY)
  INSTALL(PROGRAMS ${CMAKE_BINARY_DIR}/Install/ncanda/templates/${t} DESTINATION ${CMTK_INSTALL_DATA_DIR}/ncanda COMPONENT tools)
ENDFOREACH(t ${TEMPLATES})

