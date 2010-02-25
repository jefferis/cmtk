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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkScalarImageSimilarity.h>
#include <cmtkVolumeIO.h>

#include <math.h>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#include <iostream>

// Check ScalarImageSimilarity result against baseline.
int
testScalarImageSimilarityCheck( const char* name, const cmtk::Types::DataItem result, const cmtk::Types::DataItem baseline )
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

  int success = 0;
  
  success += testScalarImageSimilarityCheck( "GetMutualInformation", cmtk::ScalarImageSimilarity::GetMutualInformation( img0, img1 ), 1.55075 );
  success += testScalarImageSimilarityCheck( "GetNormalizedMutualInformation", cmtk::ScalarImageSimilarity::GetNormalizedMutualInformation( img0, img1 ), 1.29928 );
  success += testScalarImageSimilarityCheck( "GetRegionalMutualInformation", cmtk::ScalarImageSimilarity::GetRegionalMutualInformation( img0, img1 ), 16.8189 );
  success += testScalarImageSimilarityCheck( "GetMinusMeanSquaredDifference", cmtk::ScalarImageSimilarity::GetMinusMeanSquaredDifference( img0, img1 ), -9545.03 );
  success += testScalarImageSimilarityCheck( "GetCrossCorrelation", cmtk::ScalarImageSimilarity::GetCrossCorrelation( img0, img1 ), 0.964253 );
  success += testScalarImageSimilarityCheck( "GetGradientCorrelation", cmtk::ScalarImageSimilarity::GetGradientCorrelation( img0, img1 ), 1.8643 );
  success += testScalarImageSimilarityCheck( "GetGradientDifference", cmtk::ScalarImageSimilarity::GetGradientDifference( img0, img1 ), 22768.8 );
  success += testScalarImageSimilarityCheck( "GetPatternIntensity", cmtk::ScalarImageSimilarity::GetPatternIntensity( img0, img1 ), 98150.8 );
  success += testScalarImageSimilarityCheck( "GetDifferenceImageEntropy", cmtk::ScalarImageSimilarity::GetDifferenceImageEntropy( img0, img1 ), 2.89753 );
  success += testScalarImageSimilarityCheck( "GetCorrelationRatio", cmtk::ScalarImageSimilarity::GetCorrelationRatio( img0, img1 ), 0.933251 );

  return (success == 10) ? 0 : 1;
}

// test "ScalarImageSimilarityMemory" class
int
testScalarImageSimilarityMemory()
{
  cmtk::UniformVolume::SmartPtr testVolume( cmtk::VolumeIO::Read( CMTK_DATADIR "/spgr_3t.hdr" ) );
  if ( ! testVolume )
    {
    std::cerr << "SETUP ERROR: could not read test image 'spgr_3t.hdr'" << std::endl;
    return 1;
    }

  cmtk::ScalarImage::SmartPtr img0( testVolume->GetOrthoSlice( 2, 34 ) );
  cmtk::ScalarImage::SmartPtr img1( testVolume->GetOrthoSlice( 2, 35 ) );
  cmtk::ScalarImage::SmartPtr img2( testVolume->GetOrthoSlice( 2, 36 ) );
  
  int success = 0;

  cmtk::ScalarImageSimilarityMemory memory;
  success += testScalarImageSimilarityCheck( "GetMutualInformation#1", cmtk::ScalarImageSimilarity::GetMutualInformation( img0, img1, &memory ), 1.55075 );
  success += testScalarImageSimilarityCheck( "GetMutualInformation#2", cmtk::ScalarImageSimilarity::GetMutualInformation( img0, img1, &memory ), 1.55075 );
  success += testScalarImageSimilarityCheck( "GetMutualInformation#3", cmtk::ScalarImageSimilarity::GetMutualInformation( img1, img2, &memory ), 1.55387 );

  success += testScalarImageSimilarityCheck( "GetNormalizedMutualInformation#1", cmtk::ScalarImageSimilarity::GetNormalizedMutualInformation( img0, img1, &memory ), 1.29928 );
  success += testScalarImageSimilarityCheck( "GetNormalizedMutualInformation#2", cmtk::ScalarImageSimilarity::GetNormalizedMutualInformation( img0, img1, &memory ), 1.29928 );
  success += testScalarImageSimilarityCheck( "GetNormalizedMutualInformation#3", cmtk::ScalarImageSimilarity::GetNormalizedMutualInformation( img1, img2, &memory ), 1.30066 );

  return (success == 6) ? 0 : 1;
}
