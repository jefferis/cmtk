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

# ==========================================
# Tests for "triplanar" gui application
# (batch mode)

SET(testDriver ${CMAKE_CURRENT_BINARY_DIR}/triplanarTestDriver.sh)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/triplanarTestDriver.sh.in ${CMAKE_CURRENT_BINARY_DIR}/triplanarTestDriver.sh @ONLY)

SET(testList triplanarDefault 
	     triplanarZoom300
	     triplanarZoom25
	     triplanarSetPixelWindowLevel 
	     triplanarColormaps
	     triplanarPhantomAx
	     triplanarPhantomSa
	     triplanarPhantomCo)

IF(IGS_BUILD_NRRD OR IGS_USE_TEEM)
	SET(testList ${testList}
		     triplanarPhantomAxNrrd
		     triplanarPhantomSaNrrd
		     triplanarPhantomCoNrrd)
ENDIF(IGS_BUILD_NRRD OR IGS_USE_TEEM)

# ==========================================
# Set up all tests

FOREACH(testName ${testList})
  IF(CMTK_TESTING_MEMORYCHECK)
    ADD_TEST(${testName} /bin/sh ${testDriver} ${testName} ${MEMORYCHECK_COMMAND})
  ELSE(CMTK_TESTING_MEMORYCHECK)
    ADD_TEST(${testName} /bin/sh ${testDriver} ${testName})
  ENDIF(CMTK_TESTING_MEMORYCHECK)      
ENDFOREACH(testName ${testList})

