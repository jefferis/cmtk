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
##  $Revision: 3902 $
##      
##  $LastChangedDate: 2012-02-22 13:33:40 -0800 (Wed, 22 Feb 2012) $
##
##  $LastChangedBy: torstenrohlfing $
##      

FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.nii.gz")
INSTALL(FILES LICENSE ${files} DESTINATION ${CMTK_INSTALL_DATA_DIR}/sri24 COMPONENT runtime)

FILE(GLOB files_parc "${CMAKE_CURRENT_SOURCE_DIR}/parcellation/*.nii.gz")
INSTALL(FILES ${files_parc} DESTINATION ${CMTK_INSTALL_DATA_DIR}/sri24/parcellation COMPONENT runtime)

FILE(GLOB files_parc "${CMAKE_CURRENT_SOURCE_DIR}/fiber_masks/*.nii.gz")
INSTALL(FILES ${files_parc} DESTINATION ${CMTK_INSTALL_DATA_DIR}/sri24/fiber_masks COMPONENT runtime)
