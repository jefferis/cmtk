##
##  Copyright 2009-2011 SRI International
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
##  $Revision: 199 $
##      
##  $LastChangedDate: 2009-07-13 16:38:02 -0700 (Mon, 13 Jul 2009) $
##
##  $LastChangedBy: torstenrohlfing $
##      

PROJECT(MiniXML)
INCLUDE_REGULAR_EXPRESSION("^.*.h$")

SET(MXML_VERSION_MAJOR 2)
SET(MXML_VERSION_MINOR 6)
SET(MXML_VERSION "Mini-XML v2.6")

SET(mxml_SRCS mxml-attr.c  
	      mxml-entity.c
	      mxml-index.c
	      mxml-private.c
	      mxml-get.c
	      mxml-set.c
	      mxml-file.c
	      mxml-node.c
	      mxml-search.c
	      mxml-string.c)

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMake/CheckCInline.cmake)

CHECK_TYPE_SIZE("long long" LONG_LONG_SIZE)
IF(LONG_LONG_SIZE)
	SET(HAVE_LONG_LONG 1)
ENDIF(LONG_LONG_SIZE)

INCLUDE(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(strdup HAVE_STRDUP)
CHECK_FUNCTION_EXISTS(snprintf HAVE_SNPRINTF)
CHECK_FUNCTION_EXISTS(vsnprintf HAVE_VSNPRINTF)
CHECK_INCLUDE_FILES(pthread.h HAVE_PTHREAD_H)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

ADD_LIBRARY(cmtkMxml ${mxml_SRCS})
TARGET_LINK_LIBRARIES(cmtkMxml)

INSTALL(TARGETS cmtkMxml
  RUNTIME DESTINATION ${CMTK_INSTALL_BIN_DIR} COMPONENT runtime
  LIBRARY DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT runtime
  ARCHIVE DESTINATION ${CMTK_INSTALL_LIB_DIR} COMPONENT libraries)

FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${files} ${CMAKE_CURRENT_BINARY_DIR}/config.h DESTINATION ${CMTK_INSTALL_INCLUDE_DIR} COMPONENT headers)

SET(MXML_LIBRARIES "cmtkMxml" CACHE INTERNAL "Libraries for the mxml library")
SET(MXML_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Include directory for the mxml library headers.")
