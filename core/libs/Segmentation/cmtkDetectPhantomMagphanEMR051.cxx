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

#include <Segmentation/cmtkSphereDetectionBipolarMatchedFilterFFT.h>

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

  // create sphere detection filter based on bipolar FFT matched filtering
  SphereDetectionBipolarMatchedFilterFFT sphereDetector( *(this->m_PhantomImage) );

  // Find 1x 60mm SNR sphere
  TypedArray::SmartPtr filterResponse( sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( 0 ), 3 /*filterMargin*/ ) );
  this->m_Landmarks[0] = this->FindSphere( *filterResponse );
  this->m_Landmarks[0] = this->RefineSphereLocation( this->m_Landmarks[0], MagphanEMR051::SphereRadius( 0 ), 1 /*label*/ );
  
  // find the two 15mm spheres near estimated position
  filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( 1 ), this->GetBipolarFilterMargin() );
  
  this->m_Landmarks[1] = this->FindSphereAtDistance( *filterResponse, this->m_Landmarks[0], 90, 10 );
  this->m_Landmarks[1] = this->RefineSphereLocation( this->m_Landmarks[1], MagphanEMR051::SphereRadius( 1 ), 2 /*label*/ );

  this->m_Landmarks[2] = this->FindSphereAtDistance( *filterResponse, this->m_Landmarks[0], 60, 10 );
  this->m_Landmarks[2] = this->RefineSphereLocation( this->m_Landmarks[2], MagphanEMR051::SphereRadius( 1 ), 3 /*label*/ );

  // now use the SNR and the two 15mm spheres to define first intermediate coordinate system
  LandmarkPairList landmarkList;
  landmarkList.push_back( LandmarkPair( "SNR", MagphanEMR051::SphereCenter( 0 ), this->m_Landmarks[0] ) );
  landmarkList.push_back( LandmarkPair( "15mm@90mm", MagphanEMR051::SphereCenter( 1 ), this->m_Landmarks[1] ) );
  landmarkList.push_back( LandmarkPair( "15mm@60mm", MagphanEMR051::SphereCenter( 2 ), this->m_Landmarks[2] ) );

  // create initial rigid transformation to find approximate 10mm landmark sphere locations
  AffineXform::SmartPtr intermediateXform = FitRigidToLandmarks( landmarkList ).GetRigidXform();

  // Find 4x 30mm CNR spheres in the right order.
  filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( 1 ), this->GetBipolarFilterMargin() );
  for ( size_t i = 3; i < 7; ++i )
    {
    const Self::SpaceVectorType candidate = intermediateXform->Apply( MagphanEMR051::SphereCenter( i ) );
    this->m_Landmarks[i] = this->FindSphereAtDistance( *filterResponse, candidate, 0, 5 );
    this->m_Landmarks[i] = this->RefineSphereLocation( this->m_Landmarks[i], MagphanEMR051::SphereRadius( i ), 1+i /*label*/ );
    }

  // Find 10mm spheres in order near projected locations
  filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( 7 ), this->GetBipolarFilterMargin() );
  for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
    {
    this->m_Landmarks[i] = intermediateXform->Apply( MagphanEMR051::SphereCenter( i ) );
    this->m_Landmarks[i] = this->FindSphereAtDistance( *filterResponse, this->m_Landmarks[i], 0, 5 );
    this->m_Landmarks[i] = this->RefineSphereLocation( this->m_Landmarks[i], MagphanEMR051::SphereRadius( i ), 1+i /*label*/ );

    char name[10];
    sprintf( name, "10mm#%03d", static_cast<int>( i+1 ) );
    landmarkList.push_back( LandmarkPair( name, MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i] ) );
    }

  landmarkList.pop_front(); // remove unreliable SNR sphere before making final fit
  // create linear, not necessarily rigid, transformation based on all detected landmarks.
  this->m_PhantomToImageTransformation = FitAffineToLandmarks( landmarkList ).GetAffineXform();
  this->m_PhantomToImageTransformation->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation

  // compute residuals with current transformation
  if ( this->ComputeLandmarkFitResiduals( *(this->m_PhantomToImageTransformation) ) > this->GetLandmarkFitResidualThreshold() )
    {  
    // try to refine outliers, which probably were not properly located.
    for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i ) // we only care for the 10mm spheres.
      {
      if ( this->m_LandmarkFitResiduals[i] > this->GetLandmarkFitResidualThreshold() )
	{
	this->m_Landmarks[i] = intermediateXform->Apply( MagphanEMR051::SphereCenter( i ) );
	this->m_Landmarks[i] = this->FindSphereAtDistance( *filterResponse, this->m_Landmarks[i], 0, 0.5 * this->GetLandmarkFitResidualThreshold() );
	this->m_Landmarks[i] = this->RefineSphereLocation( this->m_Landmarks[i], MagphanEMR051::SphereRadius( i ), 1+i /*label*/ );
	}
      }
    
    LandmarkPairList landmarkList;
    landmarkList.push_back( LandmarkPair( "15mm@90mm", MagphanEMR051::SphereCenter( 1 ), this->m_Landmarks[1] ) );
    landmarkList.push_back( LandmarkPair( "15mm@60mm", MagphanEMR051::SphereCenter( 2 ), this->m_Landmarks[2] ) );
    
    for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
      {
      char name[10];
      sprintf( name, "10mm#%03d", static_cast<int>( i+1 ) );
      landmarkList.push_back( LandmarkPair( name, MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i] ) );
      }
    
    this->m_PhantomToImageTransformation = FitAffineToLandmarks( landmarkList ).GetAffineXform();
    this->m_PhantomToImageTransformation->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation
    }

  Types::Coordinate averageFittingError = 0;
  Types::Coordinate maximumFittingError = 0;
  size_t maxErrorLabel = 0;
  for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
     {
     averageFittingError += this->m_LandmarkFitResiduals[i];
     if ( this->m_LandmarkFitResiduals[i] > maximumFittingError )
       {
       maximumFittingError = this->m_LandmarkFitResiduals[i];
       maxErrorLabel = i+1;
       }
     }
  averageFittingError /= (MagphanEMR051::NumberOfSpheres-7);

  DebugOutput( 5 ) << "INFO: landmark fitting error average = " << averageFittingError << ", maximum = " <<  maximumFittingError << " maxErrLabel = " << maxErrorLabel << "\n";

  return this->m_Landmarks;
}

cmtk::Types::Coordinate
cmtk::DetectPhantomMagphanEMR051::ComputeLandmarkFitResiduals( const AffineXform& xform )
{
  Types::Coordinate maxResidual = 0;

  this->m_LandmarkFitResiduals.resize( MagphanEMR051::NumberOfSpheres );
  for ( size_t i = 0; i < MagphanEMR051::NumberOfSpheres; ++i )
    {
    this->m_LandmarkFitResiduals[i] = (this->m_Landmarks[i] - xform.Apply( MagphanEMR051::SphereCenter( i ) ) ).RootSumOfSquares();
    if ( i > 6 )
      {
      maxResidual = std::max( this->m_LandmarkFitResiduals[i], maxResidual );
      }
    }

  return maxResidual;
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

  int maxIndex = -1;
  Types::DataItem maxValue = 0;
  
  for ( size_t px = 0; px < filterResponse.GetDataSize(); ++px )
    {
    if ( (this->m_ExcludeMask->GetDataAt( px ) == 0) && (this->m_IncludeMask->GetDataAt( px ) != 0) )
      {
      const Types::DataItem value = filterResponse.ValueAt( px );
      if ( (maxIndex < 0) || (value > maxValue) )
	{
	maxValue = value;
	maxIndex = px;
	}
      }
    }

  return this->m_PhantomImage->GetGridLocation( maxIndex );
}

cmtk::DetectPhantomMagphanEMR051::SpaceVectorType 
cmtk::DetectPhantomMagphanEMR051::RefineSphereLocation( const Self::SpaceVectorType& estimate, const Types::Coordinate sphereRadius, const int label )
{
  const int margin = this->GetRefineMarginPixels();

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
  maskPainter.DrawSphere( refined, sphereRadius+this->GetSphereExcludeSafetyMargin(), label );

  return refined;
}
