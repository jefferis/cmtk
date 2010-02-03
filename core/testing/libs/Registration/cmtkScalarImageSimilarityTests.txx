/*
//
//  Copyright 2004-2010 SRI International
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
//  $Revision: 1190 $
//
//  $LastChangedDate: 2010-01-26 10:02:55 -0800 (Tue, 26 Jan 2010) $
//
//  $LastChangedBy: torsten_at_home $
//
*/

#include <cmtkScalarImageSimilarity.h>
#include <cmtkVolumeIO.h>

#include <math.h>
#include <iostream>

// Check ScalarImageSimilarity result against baseline.
bool
testScalarImageSimilarityCheck( const char* name, const cmtk::Types::DataItem result, const cmtk::Types::DataItem baseline )
{
  const cmtk::Types::DataItem tolerance = 1e-4;

  const cmtk::Types::DataItem error = 2 * fabs(result-baseline) / (fabs(result)+fabs(baseline));
  if ( !finite( result ) || (error > tolerance) )
    {
    std::cerr << name << " returned " << result << ", which exceeds tolerance for baseline " << baseline << std::endl;
    return false;
    }
  return true;
}

// test "ScalarImageSimilarity" class
int
testScalarImageSimilarity()
{
  cmtk::UniformVolume::SmartPtr testVolume( cmtk::VolumeIO::Read( CMTK_DATADIR "/spgr_3t.hdr" ) );
  if ( ! testVolume )
    {
    std::cerr << "SETUP ERROR: could not read test image 'spgr_3t.hdr'" << std::endl;
    return 1;
    }

  cmtk::ScalarImage::SmartPtr img0( testVolume->GetOrthoSlice( 2, 34 ) );
  cmtk::ScalarImage::SmartPtr img1( testVolume->GetOrthoSlice( 2, 35 ) );

  bool success = true;
  
  success = success && testScalarImageSimilarityCheck( "GetMutualInformation", cmtk::ScalarImageSimilarity::GetMutualInformation( img0, img1 ), 1.55075 );
  success = success && testScalarImageSimilarityCheck( "GetNormalizedMutualInformation", cmtk::ScalarImageSimilarity::GetNormalizedMutualInformation( img0, img1 ), 1.29928 );
  success = success && testScalarImageSimilarityCheck( "GetRegionalMutualInformation", cmtk::ScalarImageSimilarity::GetRegionalMutualInformation( img0, img1 ), 0 );
  success = success && testScalarImageSimilarityCheck( "GetMeanSquaredDifference", cmtk::ScalarImageSimilarity::GetMeanSquaredDifference( img0, img1 ), -9545.03 );
  success = success && testScalarImageSimilarityCheck( "GetCrossCorrelation", cmtk::ScalarImageSimilarity::GetCrossCorrelation( img0, img1 ), 0.964253 );
  success = success && testScalarImageSimilarityCheck( "GetGradientCorrelation", cmtk::ScalarImageSimilarity::GetGradientCorrelation( img0, img1 ), 1.8643 );
  success = success && testScalarImageSimilarityCheck( "GetGradientDifference", cmtk::ScalarImageSimilarity::GetGradientDifference( img0, img1 ), 22768.8 );
  success = success && testScalarImageSimilarityCheck( "GetPatternIntensity", cmtk::ScalarImageSimilarity::GetPatternIntensity( img0, img1 ), 98150.8 );
  success = success && testScalarImageSimilarityCheck( "GetDifferenceImageEntropy", cmtk::ScalarImageSimilarity::GetDifferenceImageEntropy( img0, img1 ), 2.89753 );
  success = success && testScalarImageSimilarityCheck( "GetCorrelationRatio", cmtk::ScalarImageSimilarity::GetCorrelationRatio( img0, img1 ), 0.933251 );

  return success ? 0 : 1;
}
