/*
//
//  Copyright 2012-2014 SRI International
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
#include <Base/cmtkUniformVolumeMorphologicalOperators.h>
#include <Base/cmtkValueSequence.h>

#include <System/cmtkDebugOutput.h>
#include <System/cmtkExitException.h>

#include <IO/cmtkVolumeIO.h>

#include <Segmentation/cmtkSphereDetectionNormalizedBipolarMatchedFilterFFT.h>
#include <Segmentation/cmtkLeastSquaresPolynomialIntensityBiasField.h>

#include <map>

cmtk::DetectPhantomMagphanEMR051::DetectPhantomMagphanEMR051( UniformVolume::SmartConstPtr& phantomImage, Self::Parameters& parameters )
  : m_Parameters( parameters ),
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
  TypedArray::SmartPtr filterResponse( sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( 0 ), this->m_Parameters.m_BipolarFilterMargin ) );
  this->m_Landmarks[0] = this->RefineSphereLocation( this->FindSphere( *filterResponse ), MagphanEMR051::SphereRadius( 0 ), 1 /*label*/ );

  // Nothing goes without the SNR sphere!
  if ( ! this->m_Landmarks[0].m_Valid )
    {
    StdErr << "ERROR: cannot find SNR sphere.\n";
    throw ExitException( 1 );
    }

  // assume that SNR sphere defines center of phantom (may not be true if phantom is broken)
  Self::SpaceVectorType phantomCenter = this->m_Landmarks[0].m_Location;
  
  // The first pass at the CNR spheres is only temporary, so store exclusion mask (SNR sphere)
  TypedArray::SmartPtr saveExcludeMaskData = this->m_ExcludeMask->GetData()->Clone();

  // Find 4x 30mm CNR spheres in the ANY order - sort them out by comparing signal intensity
  std::multimap<Types::DataItem,Self::SpaceVectorType> cnrSpheres;
  Self::SpaceVectorType cnrCenter( 0.0 );
  for ( size_t i = 3; i < 7; ++i )
    {
    (DebugOutput( 5 ) << MagphanEMR051::SphereName( i ) << "          \r").flush();
    filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( i ), this->m_Parameters.m_BipolarFilterMargin );
    const Types::Coordinate distanceFromCenter = (MagphanEMR051::SphereCenter( i ) - MagphanEMR051::SphereCenter( 0 )).RootSumOfSquares(); // at what distance from phantom center do we expect to find this sphere?
    const Self::SpaceVectorType location = this->RefineSphereLocation( this->FindSphereAtDistance( *filterResponse, this->m_Landmarks[0].m_Location, distanceFromCenter , 20 /*search band width*/ ), MagphanEMR051::SphereRadius( i ), i+1 /*label*/ );

    cnrCenter += location;
    
    Types::DataItem mean, stdev;
    this->GetSphereMeanStdDeviation( mean, stdev, location, MagphanEMR051::SphereRadius( i ), this->m_Parameters.m_ErodeCNR, 2 /*biasFieldDegree*/ );
    
    cnrSpheres.insert( std::pair<Types::DataItem,Self::SpaceVectorType>( -mean, location ) );
    }
  cnrCenter /= cnrSpheres.size();
  // the CNR spheres are only temporaty, so bring back previous exclusion mask (SNR sphere)
  this->m_ExcludeMask->SetData( saveExcludeMaskData );

  // find the two 15mm spheres near estimated position
  for ( size_t i = 1; i < 3; ++i )
    {
    (DebugOutput( 5 ) << MagphanEMR051::SphereName( i ) << "          \r").flush();
    filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( i ), this->m_Parameters.m_BipolarFilterMargin );

    const Types::Coordinate distanceFromCenter = (MagphanEMR051::SphereCenter( i ) - MagphanEMR051::SphereCenter( 0 )).RootSumOfSquares(); // at what distance from phantom center do we expect to find this sphere?
    try 
      {
      this->m_Landmarks[i] = Self::LandmarkType( this->FindSphereAtDistance( *filterResponse, this->m_Landmarks[0].m_Location, distanceFromCenter , 10 /*search band width*/ ) );
      this->m_Landmarks[i] = Self::LandmarkType( this->RefineSphereLocation( this->m_Landmarks[i].m_Location, MagphanEMR051::SphereRadius( i ), 1+i /*label*/ ) );
      }
    catch (...)
      {
      if ( !this->m_Parameters.m_TolerateTruncation )
	throw;
      }
    }

  // Check whether phantom center and the 15mm spheres actually define an orthogonal coordinate system
  if ( this->m_Landmarks[1].m_Valid && this->m_Landmarks[2].m_Valid )
    {
    // Cosine of angle between SNR sphere center and the two 15mm spheres
    const Types::Coordinate angleCosineSNR = fabs( (this->m_Landmarks[0].m_Location - this->m_Landmarks[1].m_Location).GetAngleCosine( this->m_Landmarks[0].m_Location - this->m_Landmarks[2].m_Location ) );
    // Cosine of angle between CNR 4-sphere centroid and the two 15mm spheres    
    const Types::Coordinate angleCosineCNR = fabs( (cnrCenter - this->m_Landmarks[1].m_Location).GetAngleCosine( cnrCenter - this->m_Landmarks[2].m_Location ) );

    // Use the phantom center estimate with the smallest angle cosine, i.e., the one generating the closest to an orthogonal coordinate system
    Types::Coordinate minAngleCosine;
    if ( angleCosineSNR <= angleCosineCNR )
      {
      minAngleCosine = angleCosineSNR;
      phantomCenter = this->m_Landmarks[0].m_Location;
      }
    else
      {
      // Set flag to mark that we could not use the SNR sphere.
      this->m_StatusFlags.m_FallbackCentroidCNR = true;
      this->m_StatusFlags.m_DistanceSNRtoCNR = (cnrCenter - this->m_Landmarks[0].m_Location).RootSumOfSquares();

      minAngleCosine = angleCosineCNR;
      phantomCenter = cnrCenter;
      }

    if ( (minAngleCosine > 0.087) || // more than about 5 degrees deviation from right angle - neither phantom center estimate gave an orthogonal system, this means, one or both of the 15mm spheres weren't found properly
	 ( this->m_Parameters.m_StandardOrientation &&
	   CrossProduct( this->m_Landmarks[0].m_Location - this->m_Landmarks[1].m_Location, this->m_Landmarks[0].m_Location - this->m_Landmarks[2].m_Location )[1] < 0 ) ) // wrong handedness (upside down phantom) for "standard orientation"
      {
      this->m_Landmarks[1].m_Valid = this->m_Landmarks[2].m_Valid = false;
      for ( size_t px = 0; px < this->m_ExcludeMask->GetNumberOfPixels(); ++px )
	{
	if ( this->m_ExcludeMask->GetDataAt( px ) > 1 )
	  this->m_ExcludeMask->SetDataAt( 0.0, px );
	}
      }
    }
  else
    {
    // Set flag to mark that we could not use the 15mm spheres
    this->m_StatusFlags.m_FallbackOrientationCNR = true;
    }

  // now use the SNR and the two 15mm spheres to define first intermediate coordinate system, assuming they were suiccessfully detected.
  LandmarkPairList landmarkList;
  landmarkList.push_back( LandmarkPair( "PhantomCenter", MagphanEMR051::SphereCenter( 0 ), phantomCenter ) );
  for ( size_t i = 1; i < 3; ++i )
    {
    if ( this->m_Landmarks[i].m_Valid )
      {
      landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( i ), MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i].m_Location ) );
      }
    }
    
  AffineXform::SmartPtr intermediateXform;
  // Check if all three initial landmarks were found successfully - if not, we need to fall back to CNR spheres for establishing initial gross orientation
  if (  landmarkList.size() >= 3 )
    {
    // create initial rigid transformation to find approximate 10mm landmark sphere locations
    try
      {
      intermediateXform = FitRigidToLandmarks( landmarkList ).GetRigidXform();
      }
    catch ( const AffineXform::MatrixType::SingularMatrixException& ex )
      {
      StdErr << "ERROR: singular matrix encountered while fitting intermediate rigid transformation\n";
      throw ExitException( 1 );
      }

    // Find 4x 30mm CNR spheres in the right order.
    for ( size_t i = 3; i < 7; ++i )
      {
      (DebugOutput( 5 ) << MagphanEMR051::SphereName( i ) << "          \r").flush();
      filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( i ), this->m_Parameters.m_BipolarFilterMargin );
      const Self::SpaceVectorType candidate = intermediateXform->Apply( MagphanEMR051::SphereCenter( i ) );
      this->m_Landmarks[i] = this->RefineSphereLocation( this->FindSphereAtDistance( *filterResponse, candidate, 0 /*search distance*/, 15 /*search band width*/ ), MagphanEMR051::SphereRadius( i ), 1+i /*label*/ );
      }
    }
  else
    {
    // brightest sphere is landmark #3
    this->m_Landmarks[3] = cnrSpheres.begin()->second;

    // make a temporary landmark list for the initial CNR spheres and phantom center only
    LandmarkPairList temporaryLandmarkList = landmarkList;
    temporaryLandmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( 3 ), MagphanEMR051::SphereCenter( 3 ), this->m_Landmarks[3].m_Location ) );
    
    // find sphere closest to #3 - this is #6
    for ( std::multimap<Types::DataItem,Self::SpaceVectorType>::const_iterator it = cnrSpheres.begin(); it != cnrSpheres.end(); ++it )
      {
      if ( it != cnrSpheres.begin() )
	{
	if ( (this->m_Landmarks[3].m_Location - it->second).RootSumOfSquares() < 60 ) // Other two CNR spheres are at least about 120mm away, so 60mm should be a safe threshold
	  {
	  this->m_Landmarks[6] = it->second;
	  temporaryLandmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( 6 ), MagphanEMR051::SphereCenter( 6 ), this->m_Landmarks[6].m_Location ) );
	  }
	}
      }

    // create intermediate transform based on spheres so far
    try
      {
      intermediateXform = FitRigidToLandmarks( temporaryLandmarkList ).GetRigidXform();
      }
    catch ( const AffineXform::MatrixType::SingularMatrixException& ex )
      {
      StdErr << "ERROR: singular matrix encountered while fitting intermediate rigid transformation\n";
      throw ExitException( 1 );
      }

    // Re-localize all CNR spheres, this time in the right place
    for ( size_t i = 3; i < 7; ++i )
      {
      (DebugOutput( 5 ) << MagphanEMR051::SphereName( i ) << "          \r").flush();
      filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( i ), this->m_Parameters.m_BipolarFilterMargin );
      const Self::SpaceVectorType candidate = intermediateXform->Apply( MagphanEMR051::SphereCenter( i ) );
      this->m_Landmarks[i] = this->RefineSphereLocation( this->FindSphereAtDistance( *filterResponse, candidate, 0 /*search distance*/, 15 /*search band width*/ ), MagphanEMR051::SphereRadius( i ), 1+i /*label*/ );
      landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( i ), MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i].m_Location ) );
      }

    // refine initial xform using all four CNR spheres
    intermediateXform = FitRigidToLandmarks( landmarkList ).GetRigidXform();
    }
  
  // Find 10mm spheres in order near projected locations
  for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
    {
    (DebugOutput( 5 ) << MagphanEMR051::SphereName( i ) << "          \r").flush();
    filterResponse = sphereDetector.GetFilteredImageData( MagphanEMR051::SphereRadius( i ), this->m_Parameters.m_BipolarFilterMargin );
    this->m_Landmarks[i] = Self::LandmarkType( intermediateXform->Apply( MagphanEMR051::SphereCenter( i ) ) );
    try
      {
      this->m_Landmarks[i] = Self::LandmarkType( this->FindSphereAtDistance( *filterResponse, this->m_Landmarks[i].m_Location, 0 /*search distance*/, 5 /*search band width*/ ) );
      this->m_Landmarks[i] = Self::LandmarkType( this->RefineSphereLocation( this->m_Landmarks[i].m_Location, MagphanEMR051::SphereRadius( i ), 1+i /*label*/ ) );
      }
    catch (...)
      {
      this->m_Landmarks[i].m_Valid = false;
      if ( !this->m_Parameters.m_TolerateTruncation )
	throw;
      }
    
    if ( this->m_Landmarks[i].m_Valid )
      {
      landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( i ), MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i].m_Location ) );
      // optional: refine initial xform using current landmark list
      if ( this->m_Parameters.m_RefineXformEachLandmark )
	{
	intermediateXform = FitRigidToLandmarks( landmarkList ).GetRigidXform();
	}
      }
    }

  // remove unreliable SNR sphere or CNR centroid before making final fit
  landmarkList.pop_front(); 

  // create linear, not necessarily rigid, transformation based on all detected landmarks.
  this->m_PhantomToImageTransformationAffine = FitAffineToLandmarks( landmarkList ).GetAffineXform();
  this->m_PhantomToImageTransformationAffine->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation

  this->m_PhantomToImageTransformationRigid = FitRigidToLandmarks( landmarkList ).GetRigidXform();
  this->m_PhantomToImageTransformationRigid->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation

  if ( this->m_Parameters.m_RefineOutliers )
    {
    this->RefineOutlierLandmarks( *filterResponse );
    
    landmarkList.clear();
    for ( size_t i = 1; i < 3; ++i )
      {
      if ( this->m_Landmarks[i].m_Valid )
	landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( i ), MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i].m_Location ) );
      }
    
    for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
      {
      if ( this->m_Landmarks[i].m_Valid )
	{
	landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( i ), MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i].m_Location ) );
	}
      }
    
    this->m_PhantomToImageTransformationAffine = FitAffineToLandmarks( landmarkList ).GetAffineXform();
    this->m_PhantomToImageTransformationAffine->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation
    
    this->m_PhantomToImageTransformationRigid = FitRigidToLandmarks( landmarkList ).GetRigidXform();
    this->m_PhantomToImageTransformationRigid->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation
    }

  if ( this->m_Parameters.m_ExcludeOutliers )
    {
    this->ExcludeOutlierLandmarks();

    this->m_PhantomToImageTransformationAffine = FitAffineToLandmarks( landmarkList ).GetAffineXform();
    this->m_PhantomToImageTransformationAffine->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation
    
    this->m_PhantomToImageTransformationRigid = FitRigidToLandmarks( landmarkList ).GetRigidXform();
    this->m_PhantomToImageTransformationRigid->ChangeCenter( MagphanEMR051::SphereCenter( 0 ) ); // SNR sphere center as center of rotation    
    }

  this->ComputeLandmarkFitResiduals( *(this->m_PhantomToImageTransformationAffine) );

  Types::Coordinate averageFittingError = 0;
  Types::Coordinate maximumFittingError = 0;
  size_t maxErrorLabel = 0;
  size_t validSpheres = 0;
  for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
     {
     if ( this->m_Landmarks[i].m_Valid )
       {
       averageFittingError += this->m_LandmarkFitResiduals[i];
       if ( this->m_LandmarkFitResiduals[i] > maximumFittingError )
	 {
	 maximumFittingError = this->m_LandmarkFitResiduals[i];
	 maxErrorLabel = i+1;
	 }
       ++validSpheres;
       }
     }
  averageFittingError /= validSpheres;

  DebugOutput( 2 ) << "INFO: landmark fitting error average = " << averageFittingError << " maximum = " <<  maximumFittingError << " maxErrName = " << MagphanEMR051::SphereName( maxErrorLabel-1 ) << " maxErrLabel = " << maxErrorLabel << "\n";
}

void
cmtk::DetectPhantomMagphanEMR051::RefineOutlierLandmarks( const TypedArray& filterResponse )
{
  // compute residuals with current transformation
  if ( this->ComputeLandmarkFitResiduals( *(this->m_PhantomToImageTransformationAffine) ) > this->m_Parameters.m_LandmarkFitResidualThreshold )
    {  
    // try to refine outliers, which probably were not properly located.
    for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i ) // we only care about the 10mm spheres.
      {
      if ( !this->m_Landmarks[i].m_Valid || (this->m_LandmarkFitResiduals[i] > this->m_Parameters.m_LandmarkFitResidualThreshold) )
	{
	this->m_Landmarks[i].m_Valid = false;
	const Self::SpaceVectorType predicted = this->m_PhantomToImageTransformationAffine->Apply( MagphanEMR051::SphereCenter( i ) );

	// if we predict a landmark outside the image field of view, then the phantom was not imaged completely and we need to bail
	if ( this->m_PhantomImage->IsInside( predicted ) )
	  {
	  // Find actual sphere somewhere near the predicted location
	  try
	    {
	    this->m_Landmarks[i] = this->FindSphereAtDistance( filterResponse, predicted, 0, 0.5 * this->m_Parameters.m_LandmarkFitResidualThreshold );
	    
	    // Refine detection based on local center-of-mass computation
	    const Self::SpaceVectorType refined = this->RefineSphereLocation( this->m_Landmarks[i].m_Location, MagphanEMR051::SphereRadius( i ), 1+i /*label*/ );
	    
	    // if the refined landmark is outside the image field of view, then the phantom was not imaged completely and we need to bail
	    if ( ! this->m_PhantomImage->IsInside( refined ) )
	      throw Self::OutsideFieldOfView( i, refined );
	    
	    // some spheres are darker than background - only accept refinements that improve residual fit error
	    if ( (refined - predicted).RootSumOfSquares() <= (this->m_Landmarks[i].m_Location - predicted).RootSumOfSquares() )
	      this->m_Landmarks[i] = refined;
	    }
	  catch (...)
	    {
	    if ( ! this->m_Parameters.m_TolerateTruncation )
	      throw;
	    }
	  }
	}
      }
    }
}

void
cmtk::DetectPhantomMagphanEMR051::ExcludeOutlierLandmarks()
{
  // compute residuals with current transformation
  if ( this->ComputeLandmarkFitResiduals( *(this->m_PhantomToImageTransformationAffine) ) > this->m_Parameters.m_LandmarkFitResidualThreshold )
    {  
    LandmarkPairList landmarkList;
    // add two 15mm spheres
    landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( 1 ), MagphanEMR051::SphereCenter( 1 ), this->m_Landmarks[1].m_Location ) );
    landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( 2 ), MagphanEMR051::SphereCenter( 2 ), this->m_Landmarks[2].m_Location ) );
    
    for ( size_t i = 7; i < MagphanEMR051::NumberOfSpheres; ++i )
      {
      if ( this->m_Landmarks[i].m_Valid && (this->m_LandmarkFitResiduals[i] < this->m_Parameters.m_LandmarkFitResidualThreshold) )
	{
	landmarkList.push_back( LandmarkPair( MagphanEMR051::SphereName( i ), MagphanEMR051::SphereCenter( i ), this->m_Landmarks[i].m_Location ) );
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
    if ( this->m_Landmarks[i].m_Valid )
      {
      this->m_LandmarkFitResiduals[i] = (this->m_Landmarks[i].m_Location - xform.Apply( MagphanEMR051::SphereCenter( i ) ) ).RootSumOfSquares();
      if ( i > 6 )
	{
	maxResidual = std::max( this->m_LandmarkFitResiduals[i], maxResidual );
	}
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
    if ( this->m_Landmarks[i].m_Valid )
      {
      maskPainter.DrawSphere( this->m_Landmarks[i].m_Location, MagphanEMR051::SphereRadius( i ), 1+i );
      }
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
    throw( Self::NoSphereInSearchRegion() );
    }
  
  return this->m_PhantomImage->GetGridLocation( maxIndex );
}

cmtk::DetectPhantomMagphanEMR051::SpaceVectorType 
cmtk::DetectPhantomMagphanEMR051::RefineSphereLocation( const Self::SpaceVectorType& estimate, const Types::Coordinate sphereRadius, const int label )
{
  const int margin = this->m_Parameters.m_RefineMarginPixels;

  DataGrid::IndexType centerPixelIndex;
  this->m_PhantomImage->GetClosestGridPointIndex( estimate, centerPixelIndex );

  const int nSphereRadius[3] = { margin + static_cast<int>( sphereRadius / this->m_PhantomImage->m_Delta[0] ), 
				 margin + static_cast<int>( sphereRadius / this->m_PhantomImage->m_Delta[1] ), 
				 margin + static_cast<int>( sphereRadius / this->m_PhantomImage->m_Delta[2] ) };
  
  const DataGrid::RegionType region( centerPixelIndex - DataGrid::IndexType::FromPointer( nSphereRadius ), 
				     centerPixelIndex + DataGrid::IndexType::FromPointer( nSphereRadius ) + DataGrid::IndexType( 1 ) );

  // Check whether region is entirely within image FOV
  if ( ! (region.From() >= DataGrid::IndexType( 0 ) && region.To() <= this->m_PhantomImage->m_Dims) )
    {
    throw( Self::OutsideFieldOfView( label-1, estimate ) );
    }

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

  if ( this->m_Parameters.m_CorrectSphereBiasField )
    {
    regionMask->SetData( DataGridMorphologicalOperators( regionMask ).GetEroded( 1 /*erode by 1 pixel*/ ) );
    
    std::vector<bool> regionMaskVectorErode( nPixels );
    for ( size_t i = 0; i < nPixels; ++i )
      {
      regionMaskVectorErode[i] = ( regionMask->GetDataAt( i ) != 0 );
      }

    try
      {
      regionVolume->SetData( LeastSquaresPolynomialIntensityBiasField( *regionVolume, regionMaskVectorErode, 1 /* polynomial degree */ ).GetCorrectedData() );
      }
    catch ( const LeastSquaresPolynomialIntensityBiasField::EmptyMaskException& ex )
      {
      // nothing to do; regionVolume still has its data
      }
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
  maskPainter.DrawSphere( refined, sphereRadius+this->m_Parameters.m_SphereExcludeSafetyMargin, label );

  if ( ! (MathUtil::IsFinite( refined[0] ) && MathUtil::IsFinite( refined[1] ) && MathUtil::IsFinite( refined[2] ) ) )
    {
    throw( Self::OutsideFieldOfView( label-1, estimate ) );
    }

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
    if ( includeOutliers || (this->m_LandmarkFitResiduals[i] < this->m_Parameters.m_LandmarkFitResidualThreshold) )    
      list.push_back( Landmark( MagphanEMR051::SphereName( i ), this->m_Landmarks[i].m_Location ) );
    }
  
  return list;
}

cmtk::DetectedPhantomMagphanEMR051::SmartPtr 
cmtk::DetectPhantomMagphanEMR051::GetDetectedPhantom()
{
  DetectedPhantomMagphanEMR051* detected = new DetectedPhantomMagphanEMR051( *(this->m_PhantomToImageTransformationAffine) );

  detected->m_EstimatedNonLinear = Self::SpaceVectorType( 0.0 );

  size_t countSpheresNonLinear = 0;
  const AffineXform phantomToPhysical( this->m_PhantomImage->GetImageToPhysicalMatrix() );
  for ( size_t i = 0; i < MagphanEMR051::NumberOfSpheres; ++i )
    {
    if ( this->m_Landmarks[i].m_Valid )
      {
      detected->AddLandmarkPair( MagphanEMR051::SphereName( i ), phantomToPhysical.Apply( this->m_PhantomToImageTransformationRigid->Apply( MagphanEMR051::SphereCenter( i ) ) ), phantomToPhysical.Apply( this->m_Landmarks[i].m_Location ), 
				 this->m_LandmarkFitResiduals[i], (i>=7) /*only the 10mm spheres #7 and above are considered precise enough for registration*/ );
      
      if ( i >= 7 )
	{
	const Self::SpaceVectorType residual = phantomToPhysical.Apply( this->m_PhantomToImageTransformationAffine->Apply( MagphanEMR051::SphereCenter( i ) ) ) - phantomToPhysical.Apply( this->m_Landmarks[i].m_Location );
	detected->m_EstimatedNonLinear += residual.Abs();
	++countSpheresNonLinear;
	}
      }
    }

  detected->m_EstimatedNonLinear /= countSpheresNonLinear;
  detected->m_StatusFlags = this->m_StatusFlags;

  // get SNR estimate
  Types::DataItem mean, stdev;
  this->GetSphereMeanStdDeviation( mean, stdev, this->m_Landmarks[0].m_Location, MagphanEMR051::SphereRadius( 0 ), this->m_Parameters.m_ErodeSNR, 2 /*biasFieldDegree*/ );
  detected->m_EstimatedSNR = mean / stdev;
  
  // get four CNR estimates
  for ( size_t i = 3; i < 7; ++i )
    {
    // we compute CNR per CNR sphere using formula from http://www.mr-tip.com/serv1.php?type=db1&dbs=Contrast%20to%20Noise%20Ratio (plus "fabs")
    this->GetSphereMeanStdDeviation( mean, stdev, this->m_Landmarks[i].m_Location, MagphanEMR051::SphereRadius( i ), this->m_Parameters.m_ErodeCNR, 2 /*biasFieldDegree*/ );
    detected->m_EstimatedCNR[i-3] = fabs( detected->m_EstimatedSNR - mean / stdev );
    }

  return DetectedPhantomMagphanEMR051::SmartPtr( detected );
}

void
cmtk::DetectPhantomMagphanEMR051::GetSphereMeanStdDeviation( Types::DataItem& mean, Types::DataItem& stdev, const Self::SpaceVectorType& center, const Types::Coordinate radius, const Types::Coordinate erodeBy, const int biasFieldDegree )
{
  UniformVolume::SmartPtr maskVolume( this->m_PhantomImage->CloneGrid() );
  maskVolume->CreateDataArray( TYPE_BYTE );
  maskVolume->GetData()->Fill( 0 );
  
  UniformVolumePainter maskPainter( maskVolume, UniformVolumePainter::COORDINATES_ABSOLUTE );
  maskPainter.DrawSphere( center, radius, 1 );

  if ( erodeBy )
    {
    maskVolume->SetData( UniformVolumeMorphologicalOperators( maskVolume ).GetErodedByDistance( erodeBy ) );
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
    try
      {
      dataArray = LeastSquaresPolynomialIntensityBiasField( *dataVolume, regionMaskVector, biasFieldDegree ).GetCorrectedData();
      }
    catch ( const LeastSquaresPolynomialIntensityBiasField::EmptyMaskException& ex )
      {
      // dataArray should still be the original data, but doesn't hurt to assign again, just in case.
      dataArray = dataVolume->GetData();
      }
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
