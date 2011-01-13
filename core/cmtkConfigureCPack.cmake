##
##  Copyright 1997-2010 Torsten Rohlfing
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

INCLUDE(InstallRequiredSystemLibraries)
IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  SET(CPACK_GENERATOR "ZIP;NSIS")
ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CPACK_GENERATOR "TGZ;RPM")
  ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CPACK_GENERATOR "TGZ")
  ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

SET(CPACK_SET_DESTDIR "ON")

SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CMTK -- The Computational Morphometry Toolkit")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_PACKAGE_VENDOR "SRI International - Neuroscience Program")
SET(CPACK_PACKAGE_VERSION_MAJOR "${CMTK_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${CMTK_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${CMTK_VERSION_PATCH}")
SET(CPACK_PACKAGE_FILE_NAME "CMTK-${CMTK_VERSION_MAJOR}.${CMTK_VERSION_MINOR}.${CMTK_VERSION_PATCH}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "")

set(CPACK_COMPONENTS_ALL applications gui libraries headers documentation)

INCLUDE(CPack)

INSTALL(FILES ${CPACK_RESOURCE_FILE_LICENSE} ${CPACK_PACKAGE_DESCRIPTION_FILE} DESTINATION ${CMAKE_INSTALL_PREFIX}${CMTK_INSTALL_DATA_DIR}/doc/ COMPONENT documentation)

IF(BUILD_APPS)
  CPACK_ADD_COMPONENT(applications
    DISPLAY_NAME "Command Line Applications"
    GROUP runtime)
ENDIF(BUILD_APPS)

IF(BUILD_GUI)
  CPACK_ADD_COMPONENT(gui
    DISPLAY_NAME "Graphical User Interface Applications"
    GROUP runtime)
ENDIF(BUILD_GUI)

CPACK_ADD_COMPONENT(libraries
  DISPLAY_NAME "Link Libraries"
  GROUP development)

CPACK_ADD_COMPONENT(headers
  DISPLAY_NAME "C/C++ Header Files"
  GROUP development)

IF(BUILD_DOCUMENTATION)
  CPACK_ADD_COMPONENT(documentation
    DISPLAY_NAME "CMTK Documentation"
    DISABLED)
ENDIF(BUILD_DOCUMENTATION)

CPACK_ADD_COMPONENT_GROUP(development
  DISPLAY_NAME "Development Components")

CPACK_ADD_COMPONENT_GROUP(runtime
  DISPLAY_NAME "Runtime Components")

IF(BUILD_SHARED_LIBS)
  SET(CPACK_COMPONENT_APPLICATIONS_DEPENDS libraries)
  SET(CPACK_COMPONENT_GUI_DEPENDS libraries)
ENDIF(BUILD_SHARED_LIBS)

