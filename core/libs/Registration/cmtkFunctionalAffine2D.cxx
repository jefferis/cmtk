/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkFunctionalAffine2D.h>

#include <cmtkConsole.h>

#include <cmtkScalarImage.h>
#include <cmtkTypedArraySimilarity.h>
#include <cmtkScalarImageSimilarity.h>

#include <cmtkMathUtil.h>

#include <cmtkPGM.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

FunctionalAffine2D::FunctionalAffine2D
( ScalarImage::SmartPtr& refImage, ScalarImage::SmartPtr& fltImage, const IntROI2D* fltROI ) 
  : NumberDOFs( 6 ),
    SimilarityMeasure( ScalarImageSimilarity::MI ),
    HistogramEqualization( false ),
    Parameters( 8 )
{
  RefImages.push_back( refImage );
  FltImages.push_back( fltImage );

  if ( fltROI ) 
    {
    this->FltImagesROI.push_back( ScalarImage::SmartPtr( new ScalarImage( this->FltImages[0].GetPtr(), fltROI ) ) );
    } 
  else 
    {
    this->FltImagesROI.push_back( fltImage );
    }
  
  // initialize transformation
  this->Parameters[0] = (fltROI) ? fltROI->From[0] * this->FltImages[0]->GetPixelSize(0) : 0;
  this->Parameters[1] = (fltROI) ? fltROI->From[1] * this->FltImages[0]->GetPixelSize(1) : 0;

  // no rotation
  this->Parameters[2] = 0.0;

  // no scale
  this->Parameters[3] = this->Parameters[4] = 1.0;

  // no shear
  this->Parameters[5] = 0.0;

  // center is center of ROI floating image.
  this->Parameters[6] = 0.5 * this->FltImagesROI[0]->GetPixelSize( 0 ) * ( this->FltImagesROI[0]->GetDims(0)-1 );
  this->Parameters[7] = 0.5 * this->FltImagesROI[0]->GetPixelSize( 0 ) * ( this->FltImagesROI[0]->GetDims(1)-1 );
  
  this->Transformation.Compose( this->Parameters.Elements );
}

FunctionalAffine2D::FunctionalAffine2D
( std::vector<ScalarImage::SmartPtr>& refImages, 
  std::vector<ScalarImage::SmartPtr>& fltImages,
  const IntROI2D* fltROI ) 
  : NumberDOFs( 6 ),
    SimilarityMeasure( ScalarImageSimilarity::MI ),
    HistogramEqualization( false ),
    RefImages( refImages ),
    FltImages( fltImages ),
    FltImagesROI( fltImages.size() ),
    Parameters( 8 )
{
  if ( fltROI ) 
    {
    for ( size_t i = 0; i < this->FltImages.size(); ++i )
      {
      this->FltImagesROI[i] = ScalarImage::SmartPtr( new ScalarImage( this->FltImages[i].GetPtr(), fltROI ) );
      } 
    }
  else 
    {
    for ( size_t i = 0; i < this->FltImages.size(); ++i )
      {
      this->FltImagesROI[i] = this->FltImages[i];
      }
    }
  
  // initialize transformation
  this->Parameters[0] = (fltROI) ? fltROI->From[0] * this->FltImages[0]->GetPixelSize(0) : 0;
  this->Parameters[1] = (fltROI) ? fltROI->From[1] * this->FltImages[0]->GetPixelSize(1) : 0;
  
  // no rotation
  this->Parameters[2] = 0.0;
  
  // no scale
  this->Parameters[3] = this->Parameters[4] = 1.0;

  // no shear
  this->Parameters[5] = 0.0;

  // center is center of ROI floating image.
  this->Parameters[6] = 0.5 * this->FltImagesROI[0]->GetPixelSize( 0 ) * ( this->FltImagesROI[0]->GetDims(0)-1 );
  this->Parameters[7] = 0.5 * this->FltImagesROI[0]->GetPixelSize( 0 ) * ( this->FltImagesROI[0]->GetDims(1)-1 );
  
  this->Transformation.Compose( this->Parameters.Elements );
}

void
FunctionalAffine2D::SetNumberOfBins
( const size_t minBins, const size_t maxBins )
{
  ImageSimilarityMemory.SetMinNumBins( minBins );
  if ( maxBins )
    ImageSimilarityMemory.SetMaxNumBins( maxBins );
  else
    ImageSimilarityMemory.SetMaxNumBins( minBins );
}

Types::Coordinate
FunctionalAffine2D::GetParamStep 
( const size_t idx, const Types::Coordinate mmStep ) const
{
  switch ( idx ) 
    {
    default:
    case 0:
    case 1:
      // translation
      return mmStep;
    case 2: 
    {
    // rotation
    const Types::Coordinate minSize = std::min( FltImagesROI[0]->GetDims( AXIS_X ) * FltImagesROI[0]->GetPixelSize( AXIS_X ),
					    FltImagesROI[0]->GetDims( AXIS_Y ) * FltImagesROI[0]->GetPixelSize( AXIS_Y ) );
    return MathUtil::RadToDeg( atan( 2 * mmStep / minSize ) );
    }
    case 3: 
    {
    // scale x
    const Types::Coordinate size = FltImagesROI[0]->GetDims( AXIS_X ) * FltImagesROI[0]->GetPixelSize( AXIS_X );
    return 2 * mmStep / size;
    }
    case 4: 
    {
    // scale y
    const Types::Coordinate size = FltImagesROI[0]->GetDims( AXIS_Y ) * FltImagesROI[0]->GetPixelSize( AXIS_Y );
    return 2 * mmStep / size;
    }
    case 5:
    {
    // shear
    const Types::Coordinate sizeYX = 2.0 / (FltImagesROI[0]->GetDims( AXIS_X ) * FltImagesROI[0]->GetPixelSize( AXIS_X ));
    return sizeYX;
    }
    }
  
  return 0;
}

void
FunctionalAffine2D::SetParamVector( CoordinateVector& v ) 
{
  Parameters = v;
  Transformation.Compose( Parameters.Elements );
}

void
FunctionalAffine2D::GetParamVector( CoordinateVector& v ) 
{
  v = Parameters;
}

FunctionalAffine2D::ReturnType
FunctionalAffine2D::EvaluateAt( CoordinateVector& v ) 
{
  this->SetParamVector( v );
  return this->Evaluate();
}
  
FunctionalAffine2D::ReturnType
FunctionalAffine2D::Evaluate() 
{
	Self::ReturnType result = 0;
  if ( (this->FltImagesROI.size() > 1) || (this->RefImages.size() > 1) )
    {
    std::vector<ScalarImage::SmartPtr> refImageROISmart( this->RefImages.size() );
    std::vector<const ScalarImage*> refImageROI( this->RefImages.size() );
    std::vector<const ScalarImage*> fltImageROI( this->FltImagesROI.size() );

    for ( size_t i = 0; i < RefImages.size(); ++i )
      {
      refImageROISmart[i] = ScalarImage::SmartPtr( this->RefImages[i]->InterpolateFrom( this->FltImagesROI[i], &this->Transformation ) );
      refImageROI[i] = refImageROISmart[i];
      fltImageROI[i] = this->FltImagesROI[i];
      }
    result = this->GetSimilarity( refImageROI, fltImageROI );
    }
  else
    {
    ScalarImage::SmartPtr refImageROI( this->RefImages[0]->InterpolateFrom( this->FltImagesROI[0], &this->Transformation ) );
    result = this->GetSimilarity( refImageROI, this->FltImagesROI[0] );
    }

  return result;
}

FunctionalAffine2D::ReturnType
FunctionalAffine2D::GetSimilarity
( const ScalarImage* img0,  const ScalarImage* img1 ) const
{
  switch ( this->SimilarityMeasure ) 
    {
    case ScalarImageSimilarity::MI :
      return ScalarImageSimilarity::GetMutualInformation( img0, img1, &this->ImageSimilarityMemory );
    case ScalarImageSimilarity::NMI :
      return ScalarImageSimilarity::GetNormalizedMutualInformation( img0, img1, &this->ImageSimilarityMemory );
    case ScalarImageSimilarity::RMI :
      return ScalarImageSimilarity::GetRegionalMutualInformation( img0, img1 );
    case ScalarImageSimilarity::RNMI :
      return ScalarImageSimilarity::GetRegionalMutualInformation( img0, img1, true /*NMI*/ );
    case ScalarImageSimilarity::CR :
      return ScalarImageSimilarity::GetCorrelationRatio( img0, img1 );
    case ScalarImageSimilarity::CC :
      return ScalarImageSimilarity::GetCrossCorrelation( img0, img1 );
    case ScalarImageSimilarity::MSD :
      return ScalarImageSimilarity::GetMeanSquaredDifference( img0, img1 );
    case ScalarImageSimilarity::DAE :
      return ScalarImageSimilarity::GetDifferenceImageEntropy( img0, img1 );
    case ScalarImageSimilarity::GradientCorrelation :
      return ScalarImageSimilarity::GetGradientCorrelation( img0, img1 );
    case ScalarImageSimilarity::PatternIntensity :
      return ScalarImageSimilarity::GetPatternIntensity( img0, img1 );
    default:
      return 0;
    }
}

FunctionalAffine2D::ReturnType
FunctionalAffine2D::GetSimilarity
  ( std::vector<const ScalarImage*>& imgs0,  
    std::vector<const ScalarImage*>& imgs1 ) const
{
  switch ( this->SimilarityMeasure ) 
    {
    case ScalarImageSimilarity::RMI :
      return ScalarImageSimilarity::GetMutualInformation( imgs0, imgs1 );
    case ScalarImageSimilarity::RNMI :
      return ScalarImageSimilarity::GetNormalizedMutualInformation( imgs0, imgs1 ); 
    default:
    {
    assert( imgs0.size() == imgs1.size() );
    Self::ReturnType similarity = 0;
    std::vector<const ScalarImage*>::const_iterator it0, it1;
    it0 = imgs0.begin();
    it1 = imgs1.begin();
    while ( (it0 != imgs0.end()) && (it1 != imgs1.end()) )
      {
      similarity += this->GetSimilarity( *it0, *it1 );
      ++it0;
      ++it1;
      }
    return similarity;
    }
    }
}

} // namespace cmtk
