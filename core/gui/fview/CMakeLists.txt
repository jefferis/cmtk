##
##  Copyright 1997-2009 Torsten Rohlfing
##
##  Copyright 2004-2010 SRI International
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

INCLUDE(${CMAKE_SOURCE_DIR}/libs/cmtkLibInitSrcs.cmake)

QT4_WRAP_UI(FVIEW_UIS_H fviewMainWindow.ui)

SET(fview_SRCS
  cmtkFusionViewApplication.cxx
  )

SET(fview_MOC_SRCS
  cmtkFusionViewApplication.h
  )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
  )

IF(QT_WRAP_CPP)
  QT_WRAP_CPP(fview fview_SRCS ${fview_MOC_SRCS} )
ENDIF(QT_WRAP_CPP)

ADD_EXECUTABLE(fview ${CMTK_LIB_INIT_SRCS} main.cxx ${fview_SRCS} ${FVIEW_UIS_H})
TARGET_LINK_LIBRARIES(fview ${CMTK_LIBRARIES} ${QT_LIBRARIES} ${CMTK_SYSTEM_LIBRARIES})

INSTALL(TARGETS fview RUNTIME DESTINATION bin COMPONENT gui)

