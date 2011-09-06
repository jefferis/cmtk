##
##  Copyright 1997-2009 Torsten Rohlfing
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
##  $Revision$
##
##  $LastChangedDate$
##
##  $LastChangedBy$
##

QT4_WRAP_UI(BSEDIT_UIS_H bseditMainWindow.ui)

SET(bsedit_SRCS
  cmtkBinarySegmentationEditorApplication.cxx
  )

SET(bsedit_MOC_SRCS
  cmtkBinarySegmentationEditorApplication.h
  )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
  )

IF(QT_WRAP_CPP)
  QT_WRAP_CPP(bsedit bsedit_SRCS ${bsedit_MOC_SRCS} )
ENDIF(QT_WRAP_CPP)

ADD_EXECUTABLE(bsedit main.cxx ${bsedit_SRCS} ${BSEDIT_UIS_H})
TARGET_LINK_LIBRARIES(bsedit ${CMTK_LIBRARIES} ${QT_LIBRARIES} ${CMTK_SYSTEM_LIBRARIES})

INSTALL(TARGETS bsedit RUNTIME DESTINATION bin COMPONENT gui)

