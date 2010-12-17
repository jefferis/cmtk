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

# Sources of non-templated classes.
SET(cmtkQt_SRCS
  cmtkQtIcons.cxx
  cmtkQtImageOperators.cxx
  cmtkQtProgress.cxx
  cmtkQtRenderImageRGB.cxx
  cmtkQtScrollRenderView.cxx
  cmtkQtSliderEntry.cxx
  cmtkQtTriplanarViewerBatchMode.cxx
  cmtkQtTriplanarViewer.cxx
  cmtkQtTriplanarWindow.cxx
  cmtkQtTriplanarWindowViewActions.cxx
  cmtkQtWindowLevelControls.cxx
  cmtkQtWindowLevelDialog.cxx
)

SET(cmtkQt_MOC_SRCS
  cmtkQGraphicsPixmapItemEvents.h
  cmtkQtImageOperators.h
  cmtkQtRenderImageRGB.h
  cmtkQtScrollRenderView.h
  cmtkQtSliderEntry.h
  cmtkQtTriplanarViewer.h
  cmtkQtTriplanarWindow.h
  cmtkQtWindowLevelControls.h
  cmtkQtWindowLevelDialog.h
)

#-------------------------------
# Generate Moc Files
#-------------------------------
IF(QT_WRAP_CPP)
  QT_WRAP_CPP(cmtkQt cmtkQt_SRCS ${cmtkQt_MOC_SRCS} )
ENDIF(QT_WRAP_CPP)

ADD_LIBRARY(cmtkQt ${cmtkQt_SRCS})

TARGET_LINK_LIBRARIES(cmtkQt cmtkIO cmtkBase cmtkSystem cmtkNumerics)

IF(CMTK_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(cmtkQt PROPERTIES ${CMTK_LIBRARY_PROPERTIES})
ENDIF(CMTK_LIBRARY_PROPERTIES)

INSTALL(TARGETS cmtkQt
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}${CMTK_INSTALL_BIN_DIR} COMPONENT runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}${CMTK_INSTALL_LIB_DIR} COMPONENT runtime
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}${CMTK_INSTALL_LIB_DIR} COMPONENT libraries)

FILE(GLOB files_h "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB files_txx "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${files_h} ${files_txx} DESTINATION ${CMAKE_INSTALL_PREFIX}${CMTK_INSTALL_INCLUDE_DIR}/cmtk/Qt COMPONENT headers)
