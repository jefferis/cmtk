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
##  $Revision$
##
##  $LastChangedDate$
##
##  $LastChangedBy$
##

# ==========================================
# Setup binary test drivers

SET(CMTK_LIBRARIES "cmtkRegistration;cmtkIO;cmtkBase;cmtkNumerics;cmtkSystem")

SET(DRIVERS libRegistrationTests)
IF(CMTK_BUILD_SQLITE)
  SET(DRIVERS ${DRIVERS} libRegistrationTestsSQLite)
ENDIF(CMTK_BUILD_SQLITE)

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

# ==========================================
# Tests for "cmtk::XformImageDB" class
IF(CMTK_BUILD_SQLITE)
  SET(Tests
    ImageXformDBCreate 
    ImageXformDBOpen
    ImageXformDBAddImage
    )
  
  FOREACH(T ${Tests})
    ADD_TEST(${T} ${EXECUTABLE_OUTPUT_PATH}/libRegistrationTestsSQLite ${T})
  ENDFOREACH(T Tests)
ENDIF(CMTK_BUILD_SQLITE)

