# - Try to find SQLITE3
# Once done this will define
#
#  SQLITE3_FOUND - system has SQLITE3
#  SQLITE3_INCLUDE_DIR - the SQLITE3 include directory
#  SQLITE3_LIBRARIES - Link these to use SQLITE3
#  SQLITE3_DEFINITIONS - Compiler switches required for using SQLITE3
#  SQLITE_EXECUTABLE - Comamnd line tool sqlite3
#
#  Modified from http://openlibraries.org/browser/trunk/FindSQLITE3.cmake
#
#=============================================================================
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# 
# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#   nor the names of their contributors may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
#=============================================================================

if ( SQLITE3_INCLUDE_DIR AND SQLITE3_LIBRARIES )
   # in cache already
   SET(SQLITE3_FIND_QUIETLY TRUE)
endif ( SQLITE3_INCLUDE_DIR AND SQLITE3_LIBRARIES )

# use pkg-config to get the directories and then use these values
# in the FIND_PATH() and FIND_LIBRARY() calls
if( NOT WIN32 )
  INCLUDE(FindPkgConfig)

  pkg_check_modules(SQLITE3 sqlite3)

endif( NOT WIN32 )

FIND_PATH(SQLITE3_INCLUDE_DIR NAMES sqlite3.h
  PATHS
  ${_SQLITE3IncDir}
)

FIND_LIBRARY(SQLITE3_LIBRARIES NAMES sqlite3
  PATHS
  ${_SQLITE3LinkDir}
)

INCLUDE(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set SQLITE3_FOUND to TRUE if 
# all listed variables are TRUE
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SQLite3 DEFAULT_MSG SQLITE3_LIBRARIES SQLITE3_INCLUDE_DIR)

# Extract version from header if not defined otherwise
IF(NOT SQLITE3_VERSION)
  FILE(READ ${SQLITE3_INCLUDE_DIR}/sqlite3.h SQLITE_HDR)
  STRING(REGEX MATCH "#define SQLITE_VERSION[ \t]*\"[^ \"]*\"" VERSION ${SQLITE_HDR})
  STRING(REGEX REPLACE "[^\"]*\"([^\"]*)\".*" "\\1" VERSION ${VERSION})
  SET(SQLITE3_VERSION ${VERSION} CACHE INTERNAL "SQLite3 version number")
ENDIF(NOT SQLITE3_VERSION )

MARK_AS_ADVANCED(SQLITE3_INCLUDE_DIR SQLITE3_LIBRARIES SQLITE3_EXECUTABLE)

