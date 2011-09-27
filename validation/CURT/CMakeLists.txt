##
##  Copyright 2009, 2010 SRI International
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

SET(CMTK_LIBRARIES "cmtkSegmentation;cmtkRegistration;cmtkIO;cmtkBase;cmtkNumerics;cmtkSystem")

SET(APPS curt)

FIND_PATH(ANTS_BIN_DIR ANTS)
FIND_PATH(IBSR_DIR IBSR)

OPTION(USE_QT "Use Qt (GUI applications need this)" OFF)
IF(USE_QT)
  SET (QT_MT_REQUIRED true)
  SET (QT_MIN_VERSION "4.3.0")
  FIND_PACKAGE(Qt4 COMPONENTS QtCore QtGui REQUIRED)
  IF (QT_FOUND)
    MESSAGE("Qt4 Found OK (${qt_version_str})")
    INCLUDE(${QT_USE_FILE})

    SET(APPS ${APPS} curt2d)
  ELSE(QT_FOUND)
    MESSAGE(FATAL_ERROR "No Qt4")
  ENDIF(QT_FOUND)
ENDIF(USE_QT)

FOREACH(APPLICATION ${APPS})
  ADD_EXECUTABLE(${APPLICATION} ${CMTK_LIB_INIT_SRCS} ${APPLICATION}.cxx)
  TARGET_LINK_LIBRARIES(${APPLICATION} ${QT_LIBRARIES} ${CMTK_LIBRARIES} ${CMTK_SYSTEM_LIBRARIES})
  INSTALL(TARGETS ${APPLICATION} RUNTIME DESTINATION bin)
ENDFOREACH(APPLICATION ${APPS})

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/run_validation.sh.in ${CMAKE_CURRENT_BINARY_DIR}/run_validation.sh @ONLY IMMEDIATE)
