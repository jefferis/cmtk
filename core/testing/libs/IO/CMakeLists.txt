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
##  $Revision: 1227 $
##
##  $LastChangedDate: 2010-02-08 11:45:02 -0800 (Mon, 08 Feb 2010) $
##
##  $LastChangedBy: torstenrohlfing $
##

# ==========================================
# Setup binary test drivers

SET(CMTK_LIBRARIES "cmtkIO;cmtkBase;cmtkNumerics;cmtkSystem")
SET(DRIVERS libIOTests)

FOREACH(D ${DRIVERS})
  ADD_EXECUTABLE(${D} ${D}.cxx)
  TARGET_LINK_LIBRARIES(${D} ${CMTK_LIBRARIES} ${CMTK_SYSTEM_LIBRARIES})
ENDFOREACH(D ${DRIVERS})

# ==========================================
# Tests for "cmtk::PGM" class
SET(Tests ${Tests} ReadPGM)

# ==========================================
# Tests for "cmtk::SQLite" class
SET(Tests ${Tests} 
  SQLiteCreate 
  SQLiteOpen)

FOREACH(T ${Tests})
	  ADD_TEST(${T} ${EXECUTABLE_OUTPUT_PATH}/libIOTests ${T})
ENDFOREACH(T Tests)
