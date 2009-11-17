##
##  Copyright 1997-2009 Torsten Rohlfing
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

SET(INIT_SRCS ${CMAKE_SOURCE_DIR}/libs/init/initStackBacktrace.cxx)

# Sources of non-templated classes.
SET(fusion_SRCS
  cmtkQtAdjustColormap.cxx
  cmtkQtFusionAlpha.cxx
  cmtkQtFusionEdge.cxx
  cmtkQtFusionSlicer.cxx
  cmtkQtFusionWindowTemplate.cxx
  cmtkQtSeparateView.cxx
  cmtkQtSimpleFusionApp.cxx
  cmtkQtSimpleFusionAppStudyList.cxx
  cmtkQtSimpleFusionMainWindow.cxx
  cmtkQtSimpleFusionMainWindowStudy.cxx
  cmtkQtSimpleFusionMainWindowWidget.cxx
  cmtkQtStudyNamesBox.cxx
  cmtkQtStudyWidget.cxx
  cmtkQtVolumeProperties.cxx
  main.cxx
)

SET(fusion_MOC_SRCS
  cmtkQtAdjustColormap.h
  cmtkQtFusionAlpha.h
  cmtkQtFusionEdge.h
  cmtkQtFusionSlicer.h
  cmtkQtFusionWindowTemplate.h
  cmtkQtSeparateView.h
  cmtkQtSimpleFusionApp.h
  cmtkQtSimpleFusionMainWindow.h
  cmtkQtStudyNamesBox.h
  cmtkQtStudyWidget.h
  cmtkQtVolumeProperties.h
)

#-------------------------------
# Generate Cpp/H Files
#-------------------------------
#IF(QT_WRAP_UI)
#  QT_WRAP_UI(fusion fusion_UI_H fusion_UI_CPP ${fusion_UI_SRCS} )
#ENDIF(QT_WRAP_UI)

#-------------------------------
# Generate Moc Files
#-------------------------------
IF(QT_WRAP_CPP)
  QT_WRAP_CPP(fusion fusion_SRCS ${fusion_MOC_SRCS} )
ENDIF(QT_WRAP_CPP)

SET(CMAKE_INCLUDE_CURRENT_DIR ON)
ADD_EXECUTABLE(fusion ${INIT_SRCS} ${fusion_SRCS} ${fusion_UI_CPP})
TARGET_LINK_LIBRARIES(fusion ${CMTK_LIBRARIES} ${CMTK_VTK_LIBS} ${QT_LIBRARIES} ${CMTK_SYSTEM_LIBRARIES})

INSTALL(TARGETS fusion RUNTIME DESTINATION bin)
