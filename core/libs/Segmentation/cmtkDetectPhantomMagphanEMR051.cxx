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
#include <Base/cmtkFitRigidToLandmarks.h>
#include <Base/cmtkFitAffineToLandmarks.h>
#include <Base/cmtkHistogramOtsuThreshold.h>

#include <System/cmtkDebugOutput.h>

#ifdef CMTK_USE_FFTW
#  include <Segmentation/cmtkSphereDetectionMatchedFilterFFT.h>
#endif

cmtk::DetectPhantomMagphanEMR051::DetectPhantomMagphanEMR051( UniformVolume::SmartConstPtr& phantomImage ) 
  : m_PhantomImage( phantomImage ),
    m_ExcludeMask( phantomImage->CloneGrid() ),
    m_IncludeMask( phantomImage->CloneGrid() )
{
  this->m_ExcludeMask->CreateDataArray( TYPE_BYTE, true /*setToZero*/ );
  this->m_IncludeMask->CreateDataArray( TYPE_BYTE );
}

std::vector<cmtk::DetectPhantomMagphanEMR051::SpaceVectorType> 
cmtk::DetectPhantomMagphanEMR051::GetLandmarks()
{
  this->m_Landmarks.resize( MagphanEMR051::NumberOfSpheres );

#ifdef CMTK_USE_FFTW
  SphereDetectionMatchedFilterFFT sphereDetector( *(this->m_PhantomImage) );

  // Find 1x 60mm SNR sphere
  TypedArray::SmartPtr filterResponse( sphereDetector.GetFilteredImageData( MagphanEMR051::SphereTable[0].m_Diameter / 2, 3 /*filterMargin*/ ) );
  this->m_Landmarks[0] = this->FindSphere( *filterResponse );
  this->m_Landmarks[0] = this->RefineSphereLocation( this->m_Landmarks[0], MagphanEMR051::SphereTable[0].m_Diameter / 2, 2 /*margin*/, 1 /*label*/ );
  
  // find the two 15mm spheres near estimated position
  filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereTable[1].m_Diameter / 2, 3 /*filterMargin*/ );
  
  this->m_Landmarks[1] = this->FindSphereAtDistance( *filterResponse, this->m_Landmarks[0], 90, 10 );
  this->m_Landmarks[1] = this->RefineSphereLocation( this->m_Landmarks[1], MagphanEMR051::SphereTable[1].m_Diameter / 2, 5 /*margin*/, 2 /*label*/ );

  this->m_Landmarks[2] = this->FindSphereAtDistance( *filterResponse, this->m_Landmarks[0], 60, 10 );
  this->m_Landmarks[2] = this->RefineSphereLocation( this->m_Landmarks[2], MagphanEMR051::SphereTable[2].m_Diameter / 2, 5 /*margin*/, 3 /*label*/ );

  // now use the SNR and the two 15mm spheres to define first intermediate coordinate system
  LandmarkPairList landmarkList;
  landmarkList.push_back( LandmarkPair( "SNR", Self::SpaceVectorType( MagphanEMR051::SphereTable[0].m_CenterLocation ), this->m_Landmarks[0] ) );
  landmarkList.push_back( LandmarkPair( "15mm@90mm", Self::SpaceVectorType( MagphanEMR051::SphereTable[1].m_CenterLocation ), this->m_Landmarks[1] ) );
  landmarkList.push_back( LandmarkPair( "15mm@60mm", Self::SpaceVectorType( MagphanEMR051::SphereTable[2].m_CenterLocation ), this->m_Landmarks[2] ) );

  AffineXform::SmartPtr intermediateXform = FitRigidToLandmarks( landmarkList ).GetRigidXform();

  // Find 4x 30mm CNR spheres
  filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereTable[1].m_Diameter / 2, 3 /*filterMargin*/ );
  for ( size_t i = 3; i < 7; ++i )
    {
    const Self::SpaceVectorType candidate = intermediateXform->Apply( Self::SpaceVectorType( MagphanEMR051::SphereTable[i].m_CenterLocation ) );
    this->m_Landmarks[i] = this->FindSphereAtDistance( *filterResponse, candidate, 0, 5 );
    this->m_Landmarks[i] = this->RefineSphereLocation( this->m_Landmarks[i], MagphanEMR051::SphereTable[i].m_Diameter / 2, 5 /*margin*/, 1+i /*label*/ );
    }


  // Find 10mm spheres in order near projected locations
  filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereTable[7].m_Diameter / 2, 3 /*filterMargin*/ );
  for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
    {
    const Self::SpaceVectorType candidate = intermediateXform->Apply( Self::SpaceVectorType( MagphanEMR051::SphereTable[i].m_CenterLocation ) );
    this->m_Landmarks[i] = this->FindSphereAtDistance( *filterResponse, candidate, 0, 5 );
    this->m_Landmarks[i] = this->RefineSphereLocation( this->m_Landmarks[i], MagphanEMR051::SphereTable[i].m_Diameter / 2, 5 /*margin*/, 1+i /*label*/ );

    char name[10];
    sprintf( name, "10mm#%03d", static_cast<int>( i+1 ) );
    landmarkList.push_back( LandmarkPair( name, Self::SpaceVectorType( MagphanEMR051::SphereTable[i].m_CenterLocation ), this->m_Landmarks[i] ) );
    }
#endif

  landmarkList.pop_front(); // remove unreliable SNR sphere before making final fit
  this->m_PhantomToImageTransformation = FitAffineToLandmarks( landmarkList ).GetAffineXform();
  this->m_PhantomToImageTransformation->ChangeCenter( Self::SpaceVectorType( MagphanEMR051::SphereTable[0].m_CenterLocation ) ); // SNR sphere center as center of rotation

  // compute fitting residuals
  Types::Coordinate averageFittingError = 0;
  Types::Coordinate maximumFittingError = 0;
  size_t maxErrorLabel = 0;
  for ( size_t i = 5; i < MagphanEMR051::NumberOfSpheres; ++i ) // exclude unreliable SNR and CNR spheres
    {
    const Types::Coordinate error =
      (this->m_Landmarks[i] - this->m_PhantomToImageTransformation->Apply( Self::SpaceVectorType( MagphanEMR051::SphereTable[i].m_CenterLocation ) ) ).RootSumOfSquares();
    averageFittingError += error;
    if ( error > maximumFittingError )
      {
      maximumFittingError = error;
      maxErrorLabel = i+1;
      }
    }
  averageFittingError /= (MagphanEMR051::NumberOfSpheres-5);

  DebugOutput( 5 ) << "INFO: landmark fitting error average = " << averageFittingError << ", maximum = " << maximumFittingError << " maxErrLabel = " << maxErrorLabel << "\n";

  return this->m_Landmarks;
}

cmtk::UniformVolume::SmartPtr
cmtk::DetectPhantomMagphanEMR051::GetDetectedSpheresLabelMap()
{
  // draw final sphere mask
  UniformVolumePainter maskPainter( this->m_ExcludeMask, UniformVolumePainter::COORDINATES_ABSOLUTE );
  this->m_ExcludeMask->GetData()->Fill( 0 );

  for ( size_t i = 0; i < MagphanEMR051::NumberOfSpheres; ++i )
    {
    maskPainter.DrawSphere( this->m_Landmarks[i], MagphanEMR051::SphereTable[i].m_Diameter / 2, 1+i );
    }

  return this->m_ExcludeMask;
}
    

cmtk::DetectPhantomMagphanEMR051::SpaceVectorType
cmtk::DetectPhantomMagphanEMR051::FindSphere( const TypedArray& filterResponse )
{
  size_t maxIndex = 0;
  Types::DataItem maxValue = filterResponse.ValueAt( 0 );
  
  for ( size_t px = 0; px < filterResponse.GetDataSize(); ++px )
    {
    if ( this->m_ExcludeMask->GetDataAt( px ) == 0 )
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
cmtk::DetectPhantomMagphanEMR051::FindSphereAtDistance( const TypedArray& filterResponse, const Self::SpaceVectorType& bandCenter, const Types::Coordinate bandRadius, const Types::Coordinate bandWidth )
{
  UniformVolumePainter maskPainter( this->m_IncludeMask, UniformVolumePainter::COORDINATES_ABSOLUTE );
  this->m_IncludeMask->GetData()->Fill( 0.0 );
  maskPainter.DrawSphere( bandCenter, bandRadius+bandWidth, 1 );
  if ( bandRadius > bandWidth )
    {
    maskPainter.DrawSphere( bandCenter, bandRadius-bandWidth, 0 );  
    }

  size_t maxIndex = 0;
  Types::DataItem maxValue = filterResponse.ValueAt( 0 );
  
  for ( size_t px = 0; px < filterResponse.GetDataSize(); ++px )
    {
    if ( (this->m_ExcludeMask->GetDataAt( px ) == 0) && (this->m_IncludeMask->GetDataAt( px ) != 0) )
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
cmtk::DetectPhantomMagphanEMR051::RefineSphereLocation( const Self::SpaceVectorType& estimate, const Types::Coordinate sphereRadius, const int margin, const int label )
{
  DataGrid::IndexType centerPixelIndex;
  this->m_PhantomImage->GetClosestGridPointIndex( estimate, centerPixelIndex );

  const int nSphereRadius[3] = { margin + static_cast<int>( sphereRadius / this->m_PhantomImage->m_Delta[0] ), 
				 margin + static_cast<int>( sphereRadius / this->m_PhantomImage->m_Delta[1] ), 
				 margin + static_cast<int>( sphereRadius / this->m_PhantomImage->m_Delta[2] ) };
  
  const DataGrid::RegionType region( DataGrid::IndexType( centerPixelIndex ) - DataGrid::IndexType( nSphereRadius ), 
				     DataGrid::IndexType( centerPixelIndex ) + DataGrid::IndexType( nSphereRadius ) + DataGrid::IndexType( DataGrid::IndexType::Init(1) ) );
  
  UniformVolume::SmartPtr regionVolume = this->m_PhantomImage->GetCroppedVolume( region );
  TypedArray& regionData = *(regionVolume->GetData());

  // threshold background to zero (taking it out of center-of-mass computation)
  const Types::DataItem threshold = HistogramOtsuThreshold< Histogram<unsigned int> >( *(regionData.GetHistogram( 1024 )) ).Get();
  for ( size_t i = 0; i < regionData.GetDataSize(); ++i )
    {
    if ( regionData.ValueAt( i ) < threshold )
      regionData.Set( 0.0, i );
    }
  
  const Self::SpaceVectorType refined = estimate + regionVolume->GetCenterOfMass() - regionVolume->GetCenterCropRegion();

  // update exclusion mask
  UniformVolumePainter maskPainter( this->m_ExcludeMask, UniformVolumePainter::COORDINATES_ABSOLUTE );
  maskPainter.DrawSphere( refined, sphereRadius+15, label );

  return refined;
}
