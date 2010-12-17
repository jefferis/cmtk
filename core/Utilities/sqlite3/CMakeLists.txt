##
##  Copyright 2010 Torsten Rohlfing
##
##  Copyright 2010 SRI International
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

PROJECT(sqlite3)
INCLUDE_REGULAR_EXPRESSION("^.*.h$")

SET(SQLITE_VERSION_MAJOR 3)
SET(SQLITE_VERSION_MINOR 6)
SET(SQLITE_VERSION_PATCH 22)
SET(SQLITE_VERSION "sqlite-3.6.22")

SET(sqlite_SRCS sqlite3.c)

ADD_DEFINITIONS(-DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION -DSQLITE_ENABLE_FTS3 -DSQLITE_ENABLE_RTREE)

ADD_LIBRARY(sqlite3 ${sqlite_SRCS})
TARGET_LINK_LIBRARIES(sqlite3)

INSTALL(TARGETS sqlite3
  RUNTIME DESTINATION bin COMPONENT runtime
  LIBRARY DESTINATION lib COMPONENT runtime
  ARCHIVE DESTINATION lib COMPONENT libraries)

FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}${CMTK_INSTALL_INCLUDE_DIR}/sqlite3 COMPONENT headers)

ADD_EXECUTABLE(sqlite shell.c)
TARGET_LINK_LIBRARIES(sqlite sqlite3)
INSTALL(TARGETS sqlite
  RUNTIME DESTINATION bin COMPONENT applications
  LIBRARY DESTINATION bin COMPONENT applications
  ARCHIVE DESTINATION bin COMPONENT applications)

SET(SQLITE_LIBRARIES sqlite3)
