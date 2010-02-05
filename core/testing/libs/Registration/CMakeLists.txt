##
##  Copyright 2004-2010 SRI International
##  Copyright 1997-2009 Torsten Rohlfing
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
##  $Revision: 1190 $
##
##  $LastChangedDate: 2010-01-26 10:02:55 -0800 (Tue, 26 Jan 2010) $
##
##  $LastChangedBy: torsten_at_home $
##

# ==========================================
# Setup binary test drivers

SET(CMTK_LIBRARIES "cmtkRegistration;cmtkIO;cmtkBase;cmtkNumerics;cmtkSystem")
SET(DRIVERS libRegistrationTests)

FOREACH(D ${DRIVERS})
  ADD_EXECUTABLE(${D} ${D}.cxx)
  TARGET_LINK_LIBRARIES(${D} ${CMTK_LIBRARIES} ${CMTK_SYSTEM_LIBRARIES})
ENDFOREACH(D ${DRIVERS})

# ==========================================
# Tests for "cmtk::ScalarImageSimilarity" class
SET(Tests ${Tests}
  ScalarImageSimilarity
  ScalarImageSimilarityMemory)

# ==========================================
# Tests for "cmtk::TypedArraySimilarity" class
SET(Tests ${Tests} TypedArraySimilarity)

FOREACH(T ${Tests})
  ADD_TEST(${T} ${EXECUTABLE_OUTPUT_PATH}/libRegistrationTests ${T})
ENDFOREACH(T Tests)
