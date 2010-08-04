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

#include "Registration/cmtkTypedArraySimilarity.h"
#include "IO/cmtkVolumeIO.h"

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#include <iostream>
#include <math.h>

// Check TypedArraySimilarity result against baseline.
int
testTypedArraySimilarityCheck( const char* name, const cmtk::Types::DataItem result, const cmtk::Types::DataItem baseline )
{
  const cmtk::Types::DataItem tolerance = 1e-4;

  const cmtk::Types::DataItem error = 2 * fabs(result-baseline) / (fabs(result)+fabs(baseline));
  if ( !finite( result ) || (error > tolerance) )
    {
    std::cerr << name << " returned " << result << ", which exceeds tolerance for baseline " << baseline << std::endl;
    return 0;
    }
  return 1;
}

// test "TypedArraySimilarity" class
int
testTypedArraySimilarity()
{
  cmtk::UniformVolume::SmartPtr testVolume( cmtk::VolumeIO::Read( CMTK_DATADIR "/spgr_3t.hdr" ) );
  if ( ! testVolume )
    {
    std::cerr << "SETUP ERROR: could not read test image 'spgr_3t.hdr'" << std::endl;
    return 1;
    }

  cmtk::TypedArray::SmartPtr data0( testVolume->GetOrthoSlice( 2, 34 )->GetPixelData() );
  cmtk::TypedArray::SmartPtr data1( testVolume->GetOrthoSlice( 2, 35 )->GetPixelData() );
  
  int success = 0;
  
  success += testTypedArraySimilarityCheck( "GetMutualInformation", cmtk::TypedArraySimilarity::GetMutualInformation( data0, data1 ), 1.55075 );
  success += testTypedArraySimilarityCheck( "GetNormalizedMutualInformation", cmtk::TypedArraySimilarity::GetNormalizedMutualInformation( data0, data1 ), 1.29928 );
  success += testTypedArraySimilarityCheck( "GetPeakSignalToNoiseRatio", cmtk::TypedArraySimilarity::GetPeakSignalToNoiseRatio( data0, data1 ), -4.98786 );
  success += testTypedArraySimilarityCheck( "GetMinusMeanSquaredDifference", cmtk::TypedArraySimilarity::GetMinusMeanSquaredDifference( data0, data1 ), -9545.03 );
  success += testTypedArraySimilarityCheck( "GetCrossCorrelation", cmtk::TypedArraySimilarity::GetCrossCorrelation( data0, data1 ), 0.964253 );
  success += testTypedArraySimilarityCheck( "GetCorrelationRatio", cmtk::TypedArraySimilarity::GetCorrelationRatio( data0, data1 ), 0.933251 );

  cmtk::Types::DataItem scaleFactor = 0;
  success += testTypedArraySimilarityCheck( "GetDifferenceImageEntropy", cmtk::TypedArraySimilarity::GetDifferenceArrayEntropy( data0, data1, scaleFactor ), 2.89753 );
  success += testTypedArraySimilarityCheck( "GetDifferenceImageEntropy::scaleFactor", scaleFactor, 1.01683 );

  return (success == 8) ? 0 : 1;
}
