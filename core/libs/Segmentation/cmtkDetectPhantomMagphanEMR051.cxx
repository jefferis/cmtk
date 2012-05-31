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
#include <Base/cmtkHistogramThresholdByVolume.h>
#include <Base/cmtkDataGridMorphologicalOperators.h>
#include <Base/cmtkValueSequence.h>

#include <System/cmtkDebugOutput.h>

#include <IO/cmtkVolumeIO.h>

#include <Segmentation/cmtkSphereDetectionNormalizedBipolarMatchedFilterFFT.h>
#include <Segmentation/cmtkLeastSquaresPolynomialIntensityBiasField.h>

cmtk::DetectPhantomMagphanEMR051::DetectPhantomMagphanEMR051( UniformVolume::SmartConstPtr& phantomImage ) 
  : m_CorrectSphereBiasField( true ),
    m_PhantomImage( phantomImage ),
    m_ExcludeMask( phantomImage->CloneGrid() ),
    m_IncludeMask( phantomImage->CloneGrid() )
{
  this->m_ExcludeMask->CreateDataArray( TYPE_BYTE, true /*setToZero*/ );
  this->m_IncludeMask->CreateDataArray( TYPE_BYTE );

  this->m_Landmarks.resize( MagphanEMR051::NumberOfSpheres );

  // create sphere detection filter based on bipolar FFT matched filtering
  SphereDetectionNormalizedBipolarMatchedFilterFFT sphereDetector( *(this->m_PhantomImage) );

  // Find 1x 60mm SNR sphere
  TypedArray::SmartPtr filterResponse( sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( 0 ), this->GetBipolarFilterMargin() ) );
  this->m_Landmarks[0] = this->RefineSphereLocation( this->FindSphere( *filterResponse ), MagphanEMR051::SphereRadius( 0 ), 1 /*label*/ );
  
  // find the two 15mm spheres near estimated position
  for ( size_t i = 1; i < 3; ++i )
    {
    filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( i ), this->GetBipolarFilterMargin() );

    const Types::Coordinate distanceFromCenter = (MagphanEMR051::SphereCenter( i ) - MagphanEMR051::SphereCenter( 0 )).RootSumOfSquares(); // at what distance from phantom center do we expect to find this sphere?
    this->m_Landmarks[i] = this->RefineSphereLocation( this->FindSphereAtDistance( *filterResponse, this->m_Landmarks[0], distanceFromCenter , 10 /*search band width*/ ), MagphanEMR051::SphereRadius( i ), 1+i /*label*/ );
    }

  // now use the SNR and the two 15mm spheres to define first intermediate coordinate system
  LandmarkPairList landmarkList;
  for ( size_t i = 0; i < 3; ++i )
    landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( i ), MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i] ) );

  // create initial rigid transformation to find approximate 10mm landmark sphere locations
  AffineXform::SmartPtr intermediateXform = FitRigidToLandmarks( landmarkList ).GetRigidXform();

  // Find 4x 30mm CNR spheres in the right order.
  for ( size_t i = 3; i < 7; ++i )
    {
    filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( i ), this->GetBipolarFilterMargin() );
    const Self::SpaceVectorType candidate = intermediateXform->Apply( MagphanEMR051::SphereCenter( i ) );
    this->m_Landmarks[i] = this->RefineSphereLocation( this->FindSphereAtDistance( *filterResponse, candidate, 0 /*search distance*/, 15 /*search band width*/ ), MagphanEMR051::SphereRadius( i ), 1+i /*label*/ );
    }
  
  // Find 10mm spheres in order near projected locations
  for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
    {
    filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( i ), this->GetBipolarFilterMargin() );
    this->m_Landmarks[i] = intermediateXform->Apply( MagphanEMR051::SphereCenter( i ) );
    this->m_Landmarks[i] = this->RefineSphereLocation( this->FindSphereAtDistance( *filterResponse, this->m_Landmarks[i], 0 /*search distance*/, 5 /*search band width*/ ), MagphanEMR051::SphereRadius( i ), 1+i /*label*/ );
    
    landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( i ), MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i] ) );
    }

  landmarkList.pop_front(); // remove unreliable SNR sphere before making final fit
  // create linear, not necessarily rigid, transformation based on all detected landmarks.
  this->m_PhantomToImageTransformationAffine = FitAffineToLandmarks( landmarkList ).GetAffineXform();
  this->m_PhantomToImageTransformationAffine->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation

  this->m_PhantomToImageTransformationRigid = FitRigidToLandmarks( landmarkList ).GetRigidXform();
  this->m_PhantomToImageTransformationRigid->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation

  this->RefineOutlierLandmarks( *filterResponse );
  this->ExcludeOutlierLandmarks();

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

  DebugOutput( 2 ) << "INFO: landmark fitting error average = " << averageFittingError << " maximum = " <<  maximumFittingError << " maxErrName = " << MagphanEMR051::SphereName( maxErrorLabel-1 ) << " maxErrLabel = " << maxErrorLabel << "\n";
}

void
cmtk::DetectPhantomMagphanEMR051::RefineOutlierLandmarks( const TypedArray& filterResponse )
{
  // compute residuals with current transformation
  if ( this->ComputeLandmarkFitResiduals( *(this->m_PhantomToImageTransformationAffine) ) > this->GetLandmarkFitResidualThreshold() )
    {  
    // try to refine outliers, which probably were not properly located.
    for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i ) // we only care for the 10mm spheres.
      {
      if ( this->m_LandmarkFitResiduals[i] > this->GetLandmarkFitResidualThreshold() )
	{
	const Self::SpaceVectorType predicted = this->m_Landmarks[i] = this->m_PhantomToImageTransformationAffine->Apply( MagphanEMR051::SphereCenter( i ) );
	const Self::SpaceVectorType refined = this->RefineSphereLocation( this->FindSphereAtDistance( filterResponse, this->m_Landmarks[i], 0, 0.5 * this->GetLandmarkFitResidualThreshold() ), MagphanEMR051::SphereRadius( i ), 1+i /*label*/ );
	// some spheres are darker than background - only accept refinements that improve residual fit error
	if ( (refined - predicted).RootSumOfSquares() <= (this->m_Landmarks[i] - predicted).RootSumOfSquares() )
	  this->m_Landmarks[i] = refined;
	}
      }
    
    LandmarkPairList refinedLandmarkList;
    refinedLandmarkList.push_back( LandmarkPair( "15mm@90mm", MagphanEMR051::SphereCenter( 1 ), this->m_Landmarks[1] ) );
    refinedLandmarkList.push_back( LandmarkPair( "15mm@60mm", MagphanEMR051::SphereCenter( 2 ), this->m_Landmarks[2] ) );
    
    for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
      {
      refinedLandmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( i ), MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i] ) );
      }
    }
}

void
cmtk::DetectPhantomMagphanEMR051::ExcludeOutlierLandmarks()
{
  // compute residuals with current transformation
  if ( this->ComputeLandmarkFitResiduals( *(this->m_PhantomToImageTransformationAffine) ) > this->GetLandmarkFitResidualThreshold() )
    {  
    LandmarkPairList landmarkList;
    // add two 15mm spheres
    landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( 1 ), MagphanEMR051::SphereCenter( 1 ), this->m_Landmarks[1] ) );
    landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( 2 ), MagphanEMR051::SphereCenter( 2 ), this->m_Landmarks[2] ) );
    
    for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
      {
      if ( this->m_LandmarkFitResiduals[i] < this->GetLandmarkFitResidualThreshold() )
	{
	landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( i ), MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i] ) );
	}
      }
    
    this->m_PhantomToImageTransformationAffine = FitAffineToLandmarks( landmarkList ).GetAffineXform();
    this->m_PhantomToImageTransformationAffine->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation

    this->m_PhantomToImageTransformationRigid = FitRigidToLandmarks( landmarkList ).GetRigidXform();
    this->m_PhantomToImageTransformationRigid->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation    
    }
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
    maskPainter.DrawSphere( this->m_Landmarks[i], MagphanEMR051::SphereRadius( i ), 1+i );
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
cmtk::DetectPhantomMagphanEMR051::FindSphereAtDistance
( const TypedArray& filterResponse, const Self::SpaceVectorType& bandCenter, const Types::Coordinate bandRadius, const Types::Coordinate bandWidth )
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

  if ( maxIndex < 0 )
    {
    VolumeIO::Write( *this->m_ExcludeMask, "/tmp/exclude_mask.nii" );
    throw( Self::NoSphereInSearchRegion() );
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

  UniformVolume::SmartPtr regionMask = regionVolume->CloneGrid();
  regionMask->CreateDataArray( TYPE_BYTE );
  regionMask->GetData()->Fill( 0.0 );
  
  UniformVolumePainter regionMaskPainter( regionMask, UniformVolumePainter::COORDINATES_ABSOLUTE );
  regionMaskPainter.DrawSphere( regionVolume->GetCenterCropRegion(), sphereRadius , 1 );  
  
  const size_t nPixels = regionVolume->GetNumberOfPixels();
  
  std::vector<bool> regionMaskVector( nPixels );
  for ( size_t i = 0; i < nPixels; ++i )
    {
    regionMaskVector[i] = ( regionMask->GetDataAt( i ) != 0 );
    }

  if ( this->m_CorrectSphereBiasField )
    {
    regionMask->SetData( DataGridMorphologicalOperators( regionMask ).GetEroded( 1 /*erode by 1 pixel*/ ) );
    
    std::vector<bool> regionMaskVectorErode( nPixels );
    for ( size_t i = 0; i < nPixels; ++i )
      {
      regionMaskVectorErode[i] = ( regionMask->GetDataAt( i ) != 0 );
      }
    
    regionVolume->SetData( LeastSquaresPolynomialIntensityBiasField( *regionVolume, regionMaskVectorErode, 1 /* polynomial degree */ ).GetCorrectedData() );
    }
  
  // threshold to mask out background in bias-corrected region
  const Types::DataItem threshold = HistogramOtsuThreshold< Histogram<unsigned int> >( *(regionVolume->GetData()->GetHistogram( 1024 )) ).Get();
  for ( size_t i = 0; i < nPixels; ++i )
    {
    if ( !regionMaskVector[i] || (regionVolume->GetDataAt( i ) < threshold) )
      regionVolume->SetDataAt( 0.0, i );
    }
  
  const Self::SpaceVectorType refined = estimate + regionVolume->GetCenterOfMass() - regionVolume->GetCenterCropRegion();
  
  // update exclusion mask
  UniformVolumePainter maskPainter( this->m_ExcludeMask, UniformVolumePainter::COORDINATES_ABSOLUTE );
  maskPainter.DrawSphere( refined, sphereRadius+this->GetSphereExcludeSafetyMargin(), label );

  return refined;
}

cmtk::LandmarkList
cmtk::DetectPhantomMagphanEMR051::GetExpectedLandmarks( const bool includeUnreliable ) const
{
  cmtk::LandmarkList list;

  if ( includeUnreliable )
    {
    for ( size_t i = 0; i < 7; ++i )
      {
      list.push_back( Landmark( MagphanEMR051::SphereName( i ), this->m_PhantomToImageTransformationRigid->Apply( MagphanEMR051::SphereCenter( i ) ) ) );
      }
    }
  else
    {
    // Always include 15mm spheres
    list.push_back( Landmark( MagphanEMR051::SphereName( 0 ), this->m_PhantomToImageTransformationRigid->Apply( MagphanEMR051::SphereCenter( 0 ) ) ) );
    list.push_back( Landmark( MagphanEMR051::SphereName( 1 ), this->m_PhantomToImageTransformationRigid->Apply( MagphanEMR051::SphereCenter( 1 ) ) ) );
    }
  
  for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
    {
    list.push_back( Landmark( MagphanEMR051::SphereName( i ), this->m_PhantomToImageTransformationRigid->Apply( MagphanEMR051::SphereCenter( i ) ) ) );
    }
  
  return list;
}

cmtk::LandmarkList
cmtk::DetectPhantomMagphanEMR051::GetDetectedLandmarks( const bool includeOutliers ) const
{
  cmtk::LandmarkList list;

  for ( size_t i = 0; i < MagphanEMR051::NumberOfSpheres; ++i )
    {
    if ( includeOutliers || (this->m_LandmarkFitResiduals[i] < this->GetLandmarkFitResidualThreshold()) )    
      list.push_back( Landmark( MagphanEMR051::SphereName( i ), this->m_Landmarks[i] ) );
    }
  
  return list;
}

cmtk::DetectedPhantomMagphanEMR051::SmartPtr 
cmtk::DetectPhantomMagphanEMR051::GetDetectedPhantom()
{
  DetectedPhantomMagphanEMR051* detected = new DetectedPhantomMagphanEMR051( *(this->m_PhantomToImageTransformationAffine) );

  const cmtk::AffineXform phantomToPhysical( this->m_PhantomImage->GetImageToPhysicalMatrix() );
  for ( size_t i = 0; i < MagphanEMR051::NumberOfSpheres; ++i )
    {
    detected->AddLandmarkPair( MagphanEMR051::SphereName( i ), phantomToPhysical.Apply( this->m_PhantomToImageTransformationRigid->Apply( MagphanEMR051::SphereCenter( i ) ) ), phantomToPhysical.Apply( this->m_Landmarks[i] ), 
			       this->m_LandmarkFitResiduals[i], (i>=7) /*only the 10mm spheres #7 and above are considered precise enough for registration*/ );
    }

  // get SNR estimate
  Types::DataItem mean, stdev;
  this->GetSphereMeanStdDeviation( mean, stdev, this->m_Landmarks[0], MagphanEMR051::SphereRadius( 0 ), 2 /*erodeBy*/, 2 /*biasFieldDegree*/ );
  detected->m_EstimatedSNR = mean / stdev;
  
  // get four CNR estimates
  for ( size_t i = 3; i < 7; ++i )
    {
    // we compute CNR per CNR sphere using formula from http://www.mr-tip.com/serv1.php?type=db1&dbs=Contrast%20to%20Noise%20Ratio (plus "fabs")
    this->GetSphereMeanStdDeviation( mean, stdev, this->m_Landmarks[i], MagphanEMR051::SphereRadius( i ), 2 /*erodeBy*/, 2 /*biasFieldDegree*/ );
    detected->m_EstimatedCNR[i-3] = fabs( detected->m_EstimatedSNR - mean / stdev );
    }

  return DetectedPhantomMagphanEMR051::SmartPtr( detected );
}

void
cmtk::DetectPhantomMagphanEMR051::GetSphereMeanStdDeviation( Types::DataItem& mean, Types::DataItem& stdev, const Self::SpaceVectorType& center, const Types::Coordinate radius, const int erodeBy, const int biasFieldDegree )
{
  UniformVolume::SmartPtr maskVolume( this->m_PhantomImage->CloneGrid() );
  maskVolume->CreateDataArray( TYPE_BYTE );
  maskVolume->GetData()->Fill( 0 );
  
  UniformVolumePainter maskPainter( maskVolume, UniformVolumePainter::COORDINATES_ABSOLUTE );
  maskPainter.DrawSphere( center, radius, 1 );

  if ( erodeBy )
    {
    maskVolume->SetData( DataGridMorphologicalOperators( maskVolume ).GetEroded( erodeBy ) );
    }

  // crop both mask and phantom to sphere bounding box
  UniformVolume::SmartPtr dataVolume = this->m_PhantomImage->GetCroppedVolume( maskVolume->AutoCrop( 0.5 ) );
  maskVolume = maskVolume->GetCroppedVolume();

  // make bool vector of foreground pixels
  const size_t nPixels = maskVolume->GetNumberOfPixels();
  std::vector<bool> regionMaskVector( nPixels );
  for ( size_t i = 0; i < nPixels; ++i )
    {
    regionMaskVector[i] = ( maskVolume->GetDataAt( i ) != 0 );
    }
  
  TypedArray::SmartConstPtr dataArray = dataVolume->GetData();

  // if bias correction is requested by caller, replace data with corrected data
  if ( biasFieldDegree )
    {
    dataArray = LeastSquaresPolynomialIntensityBiasField( *dataVolume, regionMaskVector, biasFieldDegree ).GetCorrectedData();
    }

  // compute summary statistics
  ValueSequence<Types::DataItem> vs;
  for ( size_t i = 0; i < nPixels; ++i )
    {
    if ( regionMaskVector[i] )
      vs.Proceed( dataArray->ValueAt( i ) );
    }

  mean = vs.GetAverage();
  stdev = sqrt( vs.GetVariance() );
}
