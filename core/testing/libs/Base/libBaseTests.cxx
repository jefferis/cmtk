/*
//
//  Copyright 2004-2010 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  This file is part of the Computational Morphometry Toolkit.
//
//  http://www.nitrc.org/projects/cmtk/
//
//  The Computational Morphometry Toolkit is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  The Computational Morphometry Toolkit is distributed in the hope that it
//  will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with the Computational Morphometry Toolkit.  If not, see
//  <http://www.gnu.org/licenses/>.
//
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <System/cmtkTestFunctionMap.h>

#include "cmtkDataGridTests.txx"
#include "cmtkEigenSystemSymmetricMatrix3x3Tests.txx"
#include "cmtkMathUtilTests.txx"
#include "cmtkMathUtilLinAlgTests.txx"
#include "cmtkParametricPlaneTests.txx"
#include "cmtkRegionTests.txx"
#include "cmtkScalarImageTests.txx"
#include "cmtkSplineWarpXformTests.txx"
#include "cmtkTypedArrayFunctionHistogramMatchingTests.txx"
#include "cmtkUniformVolumeTests.txx"

int
main( const int argc, const char* argv[] )
{
  cmtk::TestFunctionMap map;
  map.AddTest( "DataGridMatches",               &testDataGridMatches );
  map.AddTest( "EigenSystemSymmetricMatrix3x3", &testEigenSystemSymmetricMatrix3x3 );
  map.AddTest( "MathUtilEigensystem",           &testMathUtilEigensystem );
  map.AddTest( "MathUtilEigenvalues",           &testMathUtilEigenvalues );
  map.AddTest( "MathUtilUniformRandom",         &testMathUtilUniformRandom );
  map.AddTest( "ParametricPlaneMirror",         &testParametricPlaneMirror );
  map.AddTest( "ParametricPlaneMirrorOffset",   &testParametricPlaneMirrorOffset );
  map.AddTest( "RegionSizeInt",                 &testRegionSizeInt );
  map.AddTest( "RegionSizeFloat",               &testRegionSizeFloat );
  map.AddTest( "ScalarImage",                   &testScalarImage );
  map.AddTest( "SplineWarpXform",               &testSplineWarpXform );
  map.AddTest( "SplineWarpXformInverse",        &testSplineWarpXformInverse );
  map.AddTest( "TypedArrayMatchHistogram1",     &testTypedArrayMatchHistogram1 );
  map.AddTest( "TypedArrayMatchHistogram2",     &testTypedArrayMatchHistogram2 );
  map.AddTest( "TypedArrayMatchHistogram3",     &testTypedArrayMatchHistogram3 );
  map.AddTest( "TypedArrayMatchHistogram4",     &testTypedArrayMatchHistogram4 );
  map.AddTest( "UniformVolumeMatches",          &testUniformVolumeMatches );

  // is test name given on command line?
  if ( argc < 2 )
    {
    }
  else
    {
    return map.RunTestByName( argv[1] );
    }
}
