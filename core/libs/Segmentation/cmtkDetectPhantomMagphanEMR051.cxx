/*
//
//  Copyright 2012 SRI International
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

#include "cmtkDetectPhantomMagphanEMR051.h"

#include <Base/cmtkMagphanEMR051.h>
#include <Base/cmtkUniformVolumePainter.h>

#ifdef CMTK_USE_FFTW
#  include <Segmentation/cmtkSphereDetectionMatchedFilterFFT.h>
#endif

cmtk::DetectPhantomMagphanEMR051::DetectPhantomMagphanEMR051( UniformVolume::SmartConstPtr& phantomImage ) 
  : m_PhantomImage( phantomImage )
{
}

std::vector<cmtk::DetectPhantomMagphanEMR051::SpaceVectorType> 
cmtk::DetectPhantomMagphanEMR051::GetLandmarks()
{
  std::vector<Self::SpaceVectorType> landmarks( MagphanEMR051::NumberOfSpheres );

  UniformVolume::SmartPtr excludeMask( this->m_PhantomImage->CloneGrid() );
  excludeMask->CreateDataArray( TYPE_BYTE, true /*setToZero*/ );

#ifdef CMTK_USE_FFTW
  SphereDetectionMatchedFilterFFT sphereDetector( *(this->m_PhantomImage) );
#endif

  // Find 1x 60mm SNR sphere
  TypedArray::SmartPtr filterResponse( sphereDetector.GetFilteredImageData( MagphanEMR051::SphereTable[0].m_Diameter / 2, 3 /*filterMargin*/ ) );
  landmarks[0] = this->FindSphere( *filterResponse, MagphanEMR051::SphereTable[0].m_Diameter / 2, *excludeMask );
  landmarks[0] = this->RefineSphereLocation( landmarks[0], MagphanEMR051::SphereTable[0].m_Diameter / 2, 2 /*margin*/, excludeMask, 1 /*label*/ );
  
  // Find 4x 30mm CNR spheres
  filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereTable[1].m_Diameter / 2, 3 /*filterMargin*/ );
  std::vector<cmtk::DetectPhantomMagphanEMR051::SpaceVectorType> cnrLandmarks( 4 );
  for ( size_t i = 0; i < 4; ++i )
    {
    cnrLandmarks[i] = this->FindSphere( *filterResponse, MagphanEMR051::SphereTable[1+i].m_Diameter / 2, *excludeMask );
    cnrLandmarks[i] = this->RefineSphereLocation( cnrLandmarks[i], MagphanEMR051::SphereTable[1+i].m_Diameter / 2, 2 /*margin*/, excludeMask, 2+i /*label*/ );
    }

  // which of the CNR spheres is which?
  Types::DataItem averageIntensity[4] = { 0,0,0,0 };
  size_t pixelCount[4] = { 0,0,0,0 };

  for ( size_t px = 0; px < this->m_PhantomImage->GetNumberOfPixels(); ++px )
    {
    const int maskValue = static_cast<int>( excludeMask->GetDataAt( px ) );
    if ( (maskValue > 1) && (maskValue < 6) ) // there should be nothing larger than 4 right now anyway, but just to be safe...
      {
      ++pixelCount[maskValue-2];
      averageIntensity[maskValue-2] += this->m_PhantomImage->GetDataAt( px );
      }
    }
  
  for ( size_t idx = 0; idx < 4; ++idx )
    {
    averageIntensity[idx] /= pixelCount[idx];
    StdErr << averageIntensity[idx] << "\n";
    }

  for ( size_t idx = 0; idx < 4; ++idx )
    {
    size_t maxIndex = 0;
    Types::DataItem maxValue = averageIntensity[0];
    for ( size_t test = 1; test < 4; ++test )
      {
      if ( averageIntensity[test] > maxValue )
	{
	maxIndex = test;
	maxValue = averageIntensity[maxIndex];
	}
      }

    landmarks[1+idx] = cnrLandmarks[maxIndex];
    StdErr << maxIndex << "\n";
    averageIntensity[maxIndex] = 0;
    }

  // now use the SNR and the two extremal CNR spheres to define first intermediate coordinate system
  LandmarkList phantomSpaceLandmarks;
  LandmarkList imageSpaceLandmarks;

  return landmarks;
}

cmtk::DetectPhantomMagphanEMR051::SpaceVectorType
cmtk::DetectPhantomMagphanEMR051::FindSphere( const TypedArray& filterResponse, const Types::Coordinate radius, const UniformVolume& excludeMask )
{
  size_t maxIndex = 0;
  Types::DataItem maxValue = filterResponse.ValueAt( 0 );
  
  for ( size_t px = 0; px < filterResponse.GetDataSize(); ++px )
    {
    if ( excludeMask.GetDataAt( px ) == 0 )
      {
      const Types::DataItem value = filterResponse.ValueAt( px );
      if ( value > maxValue )
	{
	maxValue = value;
	maxIndex = px;
	}
      }
    }
  
  return this->m_PhantomImage->GetGridLocation( maxIndex );
}

cmtk::DetectPhantomMagphanEMR051::SpaceVectorType 
cmtk::DetectPhantomMagphanEMR051::RefineSphereLocation( const Self::SpaceVectorType& estimate, const Types::Coordinate radius, const int margin, UniformVolume::SmartPtr& excludeMask, const int label )
{
  DataGrid::IndexType centerPixelIndex;
  this->m_PhantomImage->GetClosestGridPointIndex( estimate, centerPixelIndex );

  const int nRadius[3] = { margin + static_cast<int>( radius / this->m_PhantomImage->m_Delta[0] ), 
			   margin + static_cast<int>( radius / this->m_PhantomImage->m_Delta[1] ), 
			   margin + static_cast<int>( radius / this->m_PhantomImage->m_Delta[2] ) };
  
  const DataGrid::RegionType region( DataGrid::IndexType( centerPixelIndex ) - DataGrid::IndexType( nRadius ), 
				     DataGrid::IndexType( centerPixelIndex ) + DataGrid::IndexType( nRadius ) + DataGrid::IndexType( DataGrid::IndexType::Init(1) ) );
  
  UniformVolume::SmartConstPtr regionVolume = this->m_PhantomImage->GetCroppedVolume( region );

  // threshold here
  
  const Self::SpaceVectorType refined = estimate + regionVolume->GetCenterOfMass() - regionVolume->GetCenterCropRegion();

  // update exclusion mask
  UniformVolumePainter maskPainter( excludeMask, UniformVolumePainter::COORDINATES_ABSOLUTE );
  maskPainter.DrawSphere( refined, radius, label );

  return refined;
}
