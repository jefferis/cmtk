##
##  Copyright 1997-2009 Torsten Rohlfing
##  Copyright 2004-2009 SRI International
##
##  This file is part of the Computational Morphometry Toolkit.
##
##  http:##www.nitrc.org/projects/cmtk/
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
##  <http:##www.gnu.org/licenses/>.
##
##  $Revision$
##
##  $LastChangedDate$
##
##  $LastChangedBy$
##

INCLUDE(${CMAKE_SOURCE_DIR}/libs/cmtkLibInitSrcs.cmake)

# ==========================================
# Setup binary test drivers

SET(CMTK_LIBRARIES "cmtkSystem")
SET(DRIVERS libSystemTests)

FOREACH(D ${DRIVERS})
  ADD_EXECUTABLE(${D} ${CMTK_LIB_INIT_SRCS} ${D}.cxx)
  TARGET_LINK_LIBRARIES(${D} ${CMTK_LIBRARIES})
ENDFOREACH(D ${DRIVERS})

# ==========================================
# Tests for "StackBacktrace" class
SET(Tests StackBacktrace)

FOREACH(T ${Tests})
	  ADD_TEST(${T} ${EXECUTABLE_OUTPUT_PATH}/libSystemTests ${T})
ENDFOREACH(T Tests)
