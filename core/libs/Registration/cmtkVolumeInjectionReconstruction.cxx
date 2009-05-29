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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkVolumeInjectionReconstruction.h>
#include <cmtkHistogramBase.h>
#include <cmtkVolumeIO.h>


#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

VolumeInjectionReconstruction
::VolumeInjectionReconstruction( const UniformVolume* originalImage, const int interleaveFactor, const int interleaveAxis )
  : m_NumberOfPasses( interleaveFactor ),
    m_PassWeights( interleaveFactor ),
    m_OriginalImageHistogram( new Histogram<double>( Self::NumberOfHistogramBins ) ),
    m_CorrectedImageHistogram( new Histogram<double>( Self::NumberOfHistogramBins ) )
{
  const TypedArray* originalData = originalImage->GetData();
  this->SetupHistogramKernels( originalData );

  this->m_CorrectedImage = UniformVolume::SmartPtr( originalImage->CloneGrid() );
  this->m_CorrectedImage->CreateDataArray( TYPE_FLOAT );
  
  // split original image into subvolumes and store these.
  this->m_OriginalPassImages.clear();
  for ( int pass = 0; pass < this->m_NumberOfPasses; ++pass )
    {
    UniformVolume::SmartPtr passImage( originalImage->GetInterleavedSubVolume( interleaveAxis, this->m_NumberOfPasses, pass ) );
    this->m_OriginalPassImages.push_back( passImage );
    }

  std::fill( m_PassWeights.begin(), m_PassWeights.end(), 1.0 );

  this->m_TransformationsToPassImages.clear();
  for ( int pass = 0; pass < this->m_NumberOfPasses; ++pass )
    {
    this->m_TransformationsToPassImages.push_back( AffineXform::SmartPtr( new AffineXform ) );
    }
}

VolumeInjectionReconstruction
::VolumeInjectionReconstruction( const UniformVolume* reconstructionGrid, std::vector<UniformVolume::SmartPtr>& images )
  : m_NumberOfPasses( images.size() ),
    m_PassWeights( images.size() ),
    m_OriginalImageHistogram( new Histogram<double>( Self::NumberOfHistogramBins ) ),
    m_CorrectedImageHistogram( new Histogram<double>( Self::NumberOfHistogramBins ) )
{
  const TypedArray* originalData = reconstructionGrid->GetData();
  if ( !originalData )
    originalData = images[0]->GetData();
  this->SetupHistogramKernels( originalData );

  this->m_CorrectedImage = UniformVolume::SmartPtr( reconstructionGrid->CloneGrid() );
  this->m_CorrectedImage->CreateDataArray( TYPE_FLOAT );

  this->m_OriginalPassImages = images;
  std::fill( m_PassWeights.begin(), m_PassWeights.end(), 1.0 );

  this->m_TransformationsToPassImages.clear();
  for ( int pass = 0; pass < this->m_NumberOfPasses; ++pass )
    {
    this->m_TransformationsToPassImages.push_back( AffineXform::SmartPtr( new AffineXform ) );
    }
}

int
VolumeInjectionReconstruction
::GuessInterleaveAxis
( const UniformVolume* image, const int defaultAxis )
{
  if ( (image->Dims[0] == image->Dims[1]) && (image->Dims[1] != image->Dims[2]) )
    return 2;
  if ( (image->Dims[0] == image->Dims[2]) && (image->Dims[1] != image->Dims[2]) )
    return 1;
  if ( (image->Dims[1] == image->Dims[2]) && (image->Dims[1] != image->Dims[0]) )
    return 0;

  if ( (image->Delta[0] == image->Delta[1]) && (image->Delta[1] != image->Delta[2]) )
    return 2;
  if ( (image->Delta[0] == image->Delta[2]) && (image->Delta[1] != image->Delta[2]) )
    return 1;
  if ( (image->Delta[1] == image->Delta[2]) && (image->Delta[1] != image->Delta[0]) )
    return 0;

  return defaultAxis;
}

void
VolumeInjectionReconstruction
::SetupHistogramKernels( const TypedArray* originalData )
{
  originalData->GetRange( this->m_OriginalImageMin, this->m_OriginalImageMax );
  this->m_CorrectedImageHistogram->SetRange( this->m_OriginalImageMin, this->m_OriginalImageMax );
  this->m_OriginalImageHistogram->SetRange( this->m_OriginalImageMin, this->m_OriginalImageMax );
  originalData->GetEntropy( *this->m_OriginalImageHistogram, true /*fractional*/ );
  
  const HistogramType::BinType noiseSigma = originalData->EstimateGaussianNoiseSigma( Self::NumberOfHistogramBins );
  const HistogramType::BinType kernelSigma = Self::NumberOfHistogramBins * noiseSigma / (this->m_OriginalImageMax-this->m_OriginalImageMin);
  size_t kernelRadius = static_cast<size_t>( 1 + 2 * kernelSigma );

  // We now make sure kernel radius is large enough to cover any gaps in the original image histogram.
  // This is to avoid infinite KL divergence values
  size_t runLengthZeroes = 1;
  for ( size_t i = 0; i < Self::NumberOfHistogramBins; ++i )
    {
    if ( (*this->m_OriginalImageHistogram)[i] == 0 )
      {
      ++runLengthZeroes;
      kernelRadius = std::max( kernelRadius, runLengthZeroes );
      }
    else
      {
      runLengthZeroes = 0;
      }
    }

  // Create Gaussian kernel using the previously determined sigma and cutoff radius.
  this->m_OriginalImageIntensityNoiseKernel.resize( kernelRadius );
  if ( kernelRadius > 1 )
    {
    const HistogramType::BinType normFactor = 1.0/(sqrt(2*M_PI) * kernelSigma);
    for ( size_t i = 0; i < kernelRadius; ++i )
      {
      this->m_OriginalImageIntensityNoiseKernel[i] = static_cast<HistogramType::BinType>( normFactor * exp( -MathUtil::Square( 1.0 * i / kernelSigma ) / 2 ) );
      }
    }
  else
    {
    this->m_OriginalImageIntensityNoiseKernel[0] = static_cast<HistogramType::BinType>( 1 );
    }

  // Finally, get original image histogram again, this time using the previously determined kernel.
  originalData->GetEntropy( *this->m_OriginalImageHistogram, &this->m_OriginalImageIntensityNoiseKernel[0], this->m_OriginalImageIntensityNoiseKernel.size() );
}

void
VolumeInjectionReconstruction
::ComputeTransformationsToPassImages( const int registrationMetric )
{
  this->m_TransformationsToPassImages.clear();
  
  UniformVolume::SmartPtr registerToImage = this->m_ReferenceImage ? this->m_ReferenceImage : this->m_OriginalPassImages[0];
  for ( int pass = 0; pass < this->m_NumberOfPasses; ++pass )
    {
    if ( registerToImage == this->m_OriginalPassImages[pass] )
      {
      // set identity transformation if current subvolume is used as the reference.
      this->m_TransformationsToPassImages.push_back( AffineXform::SmartPtr( new AffineXform ) );
      }
    else
      {
      // run the registration
      AffineRegistration ar;
      ar.SetVolume_1( registerToImage );
      ar.SetVolume_2( this->m_OriginalPassImages[pass] );
      
      ar.AddNumberDOFs( 6 );
      
      ar.SetInitialAlignCenters( false );
      ar.SetNoSwitch( true );
      
      ar.SetMetric( registrationMetric );
      ar.SetExploration( 4 * this->m_CorrectedImage->GetMaxDelta() );
      ar.SetAccuracy( .1 * this->m_CorrectedImage->GetMinDelta() );
      ar.SetSampling( 2 * this->m_CorrectedImage->GetMaxDelta() ); 
      
      ar.Register();
      
      this->m_TransformationsToPassImages.push_back( ar.GetTransformation()->GetInverse() );
      }
    }
}

void
VolumeInjectionReconstruction
::VolumeInjectionAnisotropic( const Types::Coordinate kernelSigmaFactor, const Types::Coordinate kernelRadiusFactor )
{
  const Types::Coordinate minusOneOverTwoSigmaSquare = -1 / (2 * kernelSigmaFactor*kernelSigmaFactor);
		    
  UniformVolume::SmartPtr& correctedImage = this->m_CorrectedImage;
  TypedArray::SmartPtr correctedImageData = correctedImage->GetData();
  const size_t correctedImageNumPixels = correctedImage->GetNumberOfPixels();

  this->m_NeighorhoodMaxPixelValues.setbounds( 1, correctedImageNumPixels );
  this->m_NeighorhoodMinPixelValues.setbounds( 1, correctedImageNumPixels );
  for ( size_t i = 1; i <= correctedImageNumPixels; ++i )
    {
    this->m_NeighorhoodMaxPixelValues(i) = this->m_OriginalImageMin;
    this->m_NeighorhoodMinPixelValues(i) = this->m_OriginalImageMax;
    }

#pragma omp parallel for schedule(dynamic)
  for ( size_t correctedPx = 0; correctedPx < correctedImageNumPixels; ++correctedPx )
    {
    double sum = 0;
    double weight = 0;

    Vector3D vCorrected;
    correctedImage->GetGridLocation( vCorrected, correctedPx );
    
    for ( int pass = 0; pass < this->m_NumberOfPasses; ++pass )
      {
      const double passImageWeight = this->m_PassWeights[pass];
      
      if ( passImageWeight > 0 )
	{
	const UniformVolume* passImage = this->m_OriginalPassImages[pass];
	const int* passImageDims = passImage->GetDims();
	const Xform* passImageXform = this->m_TransformationsToPassImages[pass];

	const Vector3D vPass = passImageXform->Apply( vCorrected );
	const Types::Coordinate passDelta[3] = { passImage->Delta[0], passImage->Delta[1], passImage->Delta[2] };

	int passGridPosition[3];
	passImage->GetVoxelIndexNoBounds( vPass, passGridPosition );

	int passGridFrom[3], passGridTo[3];
	for ( int n = 0; n < 3; ++n )
	  {
	  passGridFrom[n] = std::max( passGridPosition[n] - static_cast<int>( kernelRadiusFactor ), 0 );
	  passGridTo[n] = std::min( passGridPosition[n] + static_cast<int>( kernelRadiusFactor ) + 1, passImageDims[n] );
	  }

	Vector3D u, delta;
	for ( int k = passGridFrom[2]; k < passGridTo[2]; ++k )
	  {
	  u.XYZ[2] = passImage->GetPlaneCoord( AXIS_Z, k );
	  delta.XYZ[2] = (u.XYZ[2] - vPass.XYZ[2]) / passDelta[2];
	      
	  for ( int j = passGridFrom[1]; j < passGridTo[1]; ++j )
	    {
	    u.XYZ[1] = passImage->GetPlaneCoord( AXIS_Y, j );
	    delta.XYZ[1] = (u.XYZ[1] - vPass.XYZ[1]) / passDelta[1];
	    
	    for ( int i = passGridFrom[0]; i < passGridTo[0]; ++i )
	      {                
	      u.XYZ[0] = passImage->GetPlaneCoord( AXIS_X, i );
	      delta.XYZ[0] = (u.XYZ[0] - vPass.XYZ[0]) / passDelta[0];
	      
	      Types::DataItem passImageData;
	      if ( passImage->GetDataAt( passImageData, i, j, k ) ) 
		{
		const Types::Coordinate mahalanobis = delta.EuclidNorm();
		if ( mahalanobis <= kernelRadiusFactor )
		  {
		  const double kernelWeightPixel = passImageWeight * exp( mahalanobis*mahalanobis * minusOneOverTwoSigmaSquare );
		  
		  sum += passImageData * kernelWeightPixel;
		  weight += kernelWeightPixel;
		  
		  if ( passImageData < this->m_NeighorhoodMinPixelValues(correctedPx+1) )
		    this->m_NeighorhoodMinPixelValues(correctedPx+1) = passImageData;
		  if ( passImageData > this->m_NeighorhoodMaxPixelValues(correctedPx+1) )
		    this->m_NeighorhoodMaxPixelValues(correctedPx+1) = passImageData;
		  }
		}
	      }
	    }
	  }
	}
      }
    if ( weight > 0 )
		correctedImageData->Set( static_cast<Types::DataItem>( sum / weight ), correctedPx );
    else
      correctedImageData->SetPaddingAt( correctedPx );
    }
}

void
VolumeInjectionReconstruction
::VolumeInjectionIsotropic( const Types::Coordinate kernelSigma, const Types::Coordinate kernelRadius )
{
  const size_t correctedImageNumPixels = this->m_CorrectedImage->GetNumberOfPixels();
  const int *splattedImageDims = this->m_CorrectedImage->GetDims();

  this->m_CorrectedImage->GetData()->ClearArray();
  
  this->m_NeighorhoodMaxPixelValues.setbounds( 1, correctedImageNumPixels );
  this->m_NeighorhoodMinPixelValues.setbounds( 1, correctedImageNumPixels );
  for ( size_t i = 1; i <= correctedImageNumPixels; ++i )
    {
    this->m_NeighorhoodMaxPixelValues(i) = this->m_OriginalImageMin;
    this->m_NeighorhoodMinPixelValues(i) = this->m_OriginalImageMax;
    }

  const int kernelRadiusIndex[3] = 
  {
    1 + static_cast<int>( kernelRadius / this->m_CorrectedImage->Delta[0] ),
    1 + static_cast<int>( kernelRadius / this->m_CorrectedImage->Delta[1] ),
    1 + static_cast<int>( kernelRadius / this->m_CorrectedImage->Delta[2] )
  };

  const Types::Coordinate kernelRadiusSquare = kernelRadius * kernelRadius;
  const Types::Coordinate minusOneOverTwoSigmaSquare = -1 / (2 * kernelSigma*kernelSigma);

  std::vector<double> kernelWeights( correctedImageNumPixels );
  std::fill( kernelWeights.begin(), kernelWeights.end(), 0 );
      
  std::vector<double> splattedImage( correctedImageNumPixels );
  std::fill( splattedImage.begin(), splattedImage.end(), 0 );  
  
  for ( int pass = 0; pass < this->m_NumberOfPasses; ++pass )
    {
    const double passImageWeight = this->m_PassWeights[pass];

    if ( passImageWeight > 0 )
      {
      const UniformVolume* passImage = this->m_OriginalPassImages[pass];
      AffineXform::SmartPtr affinePassImageXform = AffineXform::SmartPtr::DynamicCastFrom( this->m_TransformationsToPassImages[pass] );
      const AffineXform* passImageXformInverse = (affinePassImageXform) ? affinePassImageXform->GetInverse() : static_cast<const AffineXform*>( NULL );
      
      const int passImageNumPixels = passImage->GetNumberOfPixels();
      for ( int offset = 0; offset < passImageNumPixels; ++offset )
	{      
	Types::DataItem passImageData;
	if ( passImage->GetDataAt( passImageData, offset ) ) 
	  {
	  int x, y, z;
	  passImage->GetIndexFromOffset( offset, x, y, z );
	  
	  Vector3D v;
	  passImage->GetGridLocation( v, x, y, z );
	  if ( passImageXformInverse )
	    {
	    passImageXformInverse->ApplyInPlace( v );
	    }
	  else
	    {
	    this->m_TransformationsToPassImages[pass]->ApplyInverseInPlace( v );
	    }
	  
	  int targetGridPosition[3];
	  if ( this->m_CorrectedImage->FindVoxel( v, targetGridPosition ) )
	    {
	    // check if neighbours are outside - if yes, set new neighbour ranges
	    int targetGridFrom[3], targetGridTo[3];
	    for ( int n = 0; n < 3; ++n )
	      {
	      targetGridFrom[n] = std::max( targetGridPosition[n] - kernelRadiusIndex[n], 0 );
	      targetGridTo[n] = std::min( targetGridPosition[n] + kernelRadiusIndex[n] + 1, splattedImageDims[n] );
	      }
	    
	    Vector3D u;
	    for ( int k = targetGridFrom[2]; k < targetGridTo[2]; ++k )
	      {
	      u.XYZ[2] = this->m_CorrectedImage->GetPlaneCoord( AXIS_Z, k );
	      
	      for ( int j = targetGridFrom[1]; j < targetGridTo[1]; ++j )
		{
		u.XYZ[1] = this->m_CorrectedImage->GetPlaneCoord( AXIS_Y, j );
		
		size_t splattedImageOffset = this->m_CorrectedImage->GetOffsetFromIndex( targetGridFrom[0], j, k );
		for ( int i = targetGridFrom[0]; i < targetGridTo[0]; ++i, ++splattedImageOffset )
		  {                
		  u.XYZ[0] = this->m_CorrectedImage->GetPlaneCoord( AXIS_X, i );
		  
		  const Types::Coordinate distanceSquare = Vector3D::SquareEuclidDistance( u, v );
		  if ( distanceSquare <= kernelRadiusSquare )
		    {
		    const double kernelWeightPixel = passImageWeight * exp( distanceSquare * minusOneOverTwoSigmaSquare );
		    
		    splattedImage[splattedImageOffset] += passImageData * kernelWeightPixel;
		    kernelWeights[splattedImageOffset] += kernelWeightPixel;
		    
		    if ( passImageData < this->m_NeighorhoodMinPixelValues(splattedImageOffset+1) )
		      this->m_NeighorhoodMinPixelValues(splattedImageOffset+1) = passImageData;
		    if ( passImageData > this->m_NeighorhoodMaxPixelValues(splattedImageOffset+1) )
		      this->m_NeighorhoodMaxPixelValues(splattedImageOffset+1) = passImageData;
		    }
		  }
		}
	      }
	    }
	  }
	}      
      }
    }
  
#pragma omp parallel for
  for ( size_t idx = 0; idx < correctedImageNumPixels; ++idx )
    {
    const Types::DataItem kernelWeightPixel = static_cast<Types::DataItem>( kernelWeights[idx] );
    if ( kernelWeightPixel > 0 ) // check if pixel is a neighbour
      {
		  this->m_CorrectedImage->SetDataAt( static_cast<Types::DataItem>( splattedImage[idx] / kernelWeightPixel ), idx ); // Set normalized data on grid
      }
    }
}

UniformVolume::SmartPtr&
VolumeInjectionReconstruction
::GetCorrectedImage()
{
  return this->m_CorrectedImage;
}

void
VolumeInjectionReconstruction
::SetReferenceImage( UniformVolume::SmartPtr& referenceImage )
{
  this->m_ReferenceImage = referenceImage;
}

double
VolumeInjectionReconstruction
::GetOriginalToCorrectedImageKLD( const ap::real_1d_array& x )
{
  this->m_CorrectedImageHistogram->Reset();
  for ( int i = x.getlowbound(); i <= x.gethighbound(); ++i )
    this->m_CorrectedImageHistogram->AddWeightedSymmetricKernel
	( this->m_CorrectedImageHistogram->ValueToBin( static_cast<Types::DataItem>( x(i) ) ), this->m_OriginalImageIntensityNoiseKernel.size(), &this->m_OriginalImageIntensityNoiseKernel[0] );
  const double kld = this->m_CorrectedImageHistogram->GetKullbackLeiblerDivergence( *this->m_OriginalImageHistogram );

  return kld;
}

double
VolumeInjectionReconstruction
::ComputeCorrectedImageLaplacianNorm( const ap::real_1d_array& correctedImagePixels )
{
  const UniformVolume* correctedImage = this->m_CorrectedImage;
  const size_t correctedImageNumPixels = correctedImage->GetNumberOfPixels();
  this->m_CorrectedImageLaplacians.resize( correctedImageNumPixels );

  const int* correctedImageDims = correctedImage->GetDims();
  const int nextI = 1;
  const int nextJ = nextI * correctedImageDims[0];
  const int nextK = nextJ * correctedImageDims[1];

  double lnorm = 0;
#pragma omp parallel for reduction(+:lnorm)
  for ( size_t idx = 1; idx <= correctedImageNumPixels; ++idx )
    {
    int x, y, z;
    correctedImage->GetIndexFromOffset( idx-1, x, y, z );
    
    const int xm = (x>0) ? idx-nextI : idx+nextI;
    const int ym = (y>0) ? idx-nextJ : idx+nextJ;
    const int zm = (z>0) ? idx-nextK : idx+nextK;
    
    const int xp = (x+1 < correctedImageDims[0]) ? idx+nextI : idx-nextI;
    const int yp = (y+1 < correctedImageDims[1]) ? idx+nextJ : idx-nextJ;
    const int zp = (z+1 < correctedImageDims[2]) ? idx+nextK : idx-nextK;
    
    const double l = 
      correctedImagePixels( xm ) + correctedImagePixels( xp ) +
      correctedImagePixels( ym ) + correctedImagePixels( yp ) +
      correctedImagePixels( zm ) + correctedImagePixels( zp ) - 
      6 * correctedImagePixels( idx );
    
    this->m_CorrectedImageLaplacians[idx-1] = l;
    lnorm += l*l;
    }
  return lnorm / correctedImageNumPixels;
}

void
VolumeInjectionReconstruction
::AddLaplacianGradientImage( ap::real_1d_array& g, const ap::real_1d_array&, const double weight ) const
{
  const UniformVolume* correctedImage = this->m_CorrectedImage;
  const size_t correctedImageNumPixels = correctedImage->GetNumberOfPixels();
  const int* correctedImageDims = correctedImage->GetDims();

  const int nextI = 1;
  const int nextJ = nextI * correctedImageDims[0];
  const int nextK = nextJ * correctedImageDims[1];
  
#pragma omp parallel for
  for ( size_t idx = 0; idx < correctedImageNumPixels; ++idx )
    {
    int x, y, z;
    correctedImage->GetIndexFromOffset( idx, x, y, z );
    
    const int xm = (x>0) ? idx-nextI : idx+nextI;
    const int ym = (y>0) ? idx-nextJ : idx+nextJ;
    const int zm = (z>0) ? idx-nextK : idx+nextK;

    const int xp = (x+1 < correctedImageDims[0]) ? idx+nextI : idx-nextI;
    const int yp = (y+1 < correctedImageDims[1]) ? idx+nextJ : idx-nextJ;
    const int zp = (z+1 < correctedImageDims[2]) ? idx+nextK : idx-nextK;

    g(idx+1) += (2 * weight / correctedImageNumPixels) *
      ( this->m_CorrectedImageLaplacians[xm] + this->m_CorrectedImageLaplacians[xp] +
	this->m_CorrectedImageLaplacians[ym] + this->m_CorrectedImageLaplacians[yp] +
	this->m_CorrectedImageLaplacians[zm] + this->m_CorrectedImageLaplacians[zp] - 
	6 * this->m_CorrectedImageLaplacians[idx] );
    }

#ifdef __IGNORE__
  UniformVolume::SmartPtr image( this->m_CorrectedImage->CloneGrid() );
  image->CreateDataArray( TYPE_FLOAT );
  
  for ( size_t i = 0; i < correctedImageNumPixels; ++i )
    image->SetDataAt( g(i+1), i );
  VolumeIO::Write( image, "gradient.nrrd" );

  for ( size_t i = 0; i < correctedImageNumPixels; ++i )
    image->SetDataAt( this->m_CorrectedImageLaplacians[i], i );
  VolumeIO::Write( image, "laplacian.nrrd" );

  VolumeIO::Write( this->m_CorrectedImage, "corrected.nrrd" );

  exit(1);
#endif
}

} // namespace cmtk
