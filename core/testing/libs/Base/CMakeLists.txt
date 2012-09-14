##
##  Copyright 1997-2009 Torsten Rohlfing
##
##  Copyright 2004-2012 SRI International
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

SET(CMTK_LIBRARIES "cmtkBase;cmtkNumerics;cmtkSystem")
SET(DRIVERS libBaseTests)

FOREACH(D ${DRIVERS})
  ADD_EXECUTABLE(${D} ${D}.cxx)
  TARGET_LINK_LIBRARIES(${D} ${CMTK_LIBRARIES} ${CMTK_SYSTEM_LIBRARIES})
ENDFOREACH(D ${DRIVERS})

# ==========================================
# Tests for "cmtk::DataGrid" class
SET(Tests ${Tests} DataGridMatches)

# ==========================================
# Tests for "cmtk::Region" class
SET(Tests ${Tests} 
  RegionSizeInt
  RegionSizeFloat)

# ==========================================
# Tests for "cmtk::ParametricPlane" class
SET(Tests ${Tests} 
  ParametricPlaneMirror
  ParametricPlaneMirrorOffset)

# ==========================================
# Tests for "cmtk::ScalarImage" class
SET(Tests ${Tests} ScalarImage)

# ==========================================
# Tests for "cmtk::SplineWarpXform" class
SET(Tests ${Tests} SplineWarpXform SplineWarpXformInverse)

# ==========================================
# Tests for "cmtk::SymmetricMatrix" class
SET(Tests ${Tests} 
  SymmetricMatrix
  SymmetricMatrixResize
  SymmetricMatrixEqual)

# ==========================================
# Tests for "cmtk::TypedArray" class
SET(Tests ${Tests} 
  TypedArrayMatchHistogram1
  TypedArrayMatchHistogram2
  TypedArrayMatchHistogram3
  TypedArrayMatchHistogram4)

# ==========================================
# Tests for "cmtk::UniformVolume" class
SET(Tests ${Tests} UniformVolumeMatches)

# ==========================================
# Tests for "cmtk::MathUtil" functions
SET(Tests ${Tests} 
	  EigenSystemSymmetricMatrix3x3
	  MathUtilUniformRandom)

FOREACH(T ${Tests})
	  ADD_TEST(${T} ${EXECUTABLE_OUTPUT_PATH}/libBaseTests ${T})
ENDFOREACH(T Tests)
