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

#include "cmtkSphereDetectionNormalizedBipolarMatchedFilterFFT.h"

cmtk::SphereDetectionNormalizedBipolarMatchedFilterFFT::SphereDetectionNormalizedBipolarMatchedFilterFFT( const UniformVolume& image )
  :   m_NumberOfPixels( image.GetNumberOfPixels() ),
      m_ImageDims( image.m_Dims ),
      m_PixelSize( image.m_Delta ),
      m_SphereRadius( 0 ),
      m_MarginWidth( -1 )
{
  this->m_ImageFT = fftw_alloc_complex( this->m_NumberOfPixels );
  this->m_ImageSquareFT = fftw_alloc_complex( this->m_NumberOfPixels );

  this->m_FilterFT = fftw_alloc_complex( this->m_NumberOfPixels );
  this->m_FilterSquareFT = fftw_alloc_complex( this->m_NumberOfPixels );
  this->m_FilterMaskFT = fftw_alloc_complex( this->m_NumberOfPixels );
  this->m_FilterMaskFT2 = fftw_alloc_complex( this->m_NumberOfPixels );

  this->m_PlanFilter = fftw_plan_dft_3d( this->m_ImageDims[2], this->m_ImageDims[1], this->m_ImageDims[0], this->m_FilterFT, this->m_FilterFT, FFTW_FORWARD, FFTW_ESTIMATE );
  this->m_PlanFilterSquare = fftw_plan_dft_3d( this->m_ImageDims[2], this->m_ImageDims[1], this->m_ImageDims[0], this->m_FilterSquareFT, this->m_FilterSquareFT, FFTW_FORWARD, FFTW_ESTIMATE );
  this->m_PlanFilterMask = fftw_plan_dft_3d( this->m_ImageDims[2], this->m_ImageDims[1], this->m_ImageDims[0], this->m_FilterMaskFT, this->m_FilterMaskFT, FFTW_FORWARD, FFTW_ESTIMATE );

  this->m_PlanBackward = fftw_plan_dft_3d( this->m_ImageDims[2], this->m_ImageDims[1], this->m_ImageDims[0], this->m_FilterFT, this->m_FilterFT, FFTW_BACKWARD, FFTW_ESTIMATE );
  this->m_PlanBackwardMask = fftw_plan_dft_3d( this->m_ImageDims[2], this->m_ImageDims[1], this->m_ImageDims[0], this->m_FilterMaskFT, this->m_FilterMaskFT, FFTW_BACKWARD, FFTW_ESTIMATE );
  this->m_PlanBackwardMask2 = fftw_plan_dft_3d( this->m_ImageDims[2], this->m_ImageDims[1], this->m_ImageDims[0], this->m_FilterMaskFT2, this->m_FilterMaskFT2, FFTW_BACKWARD, FFTW_ESTIMATE );
  
  // initialize image FT
  fftw_plan plan_image = fftw_plan_dft_3d( this->m_ImageDims[2], this->m_ImageDims[1], this->m_ImageDims[0], this->m_ImageFT, this->m_ImageFT, FFTW_FORWARD, FFTW_ESTIMATE );
  fftw_plan plan_image_square = fftw_plan_dft_3d( this->m_ImageDims[2], this->m_ImageDims[1], this->m_ImageDims[0], this->m_ImageSquareFT, this->m_ImageSquareFT, FFTW_FORWARD, FFTW_ESTIMATE );
  for ( size_t n = 0; n < this->m_NumberOfPixels; ++n )
    {
    this->m_ImageFT[n][0] = image.GetDataAt( n );
    this->m_ImageFT[n][1] = 0;

    this->m_ImageSquareFT[n][0] = this->m_ImageFT[n][0] * this->m_ImageFT[n][0];
    this->m_ImageSquareFT[n][1] = 0;
    }
  
  fftw_execute( plan_image );
  fftw_execute( plan_image_square );

  fftw_destroy_plan( plan_image );
  fftw_destroy_plan( plan_image_square );
}

cmtk::SphereDetectionNormalizedBipolarMatchedFilterFFT::~SphereDetectionNormalizedBipolarMatchedFilterFFT()
{
  fftw_destroy_plan( this->m_PlanBackward );
  fftw_destroy_plan( this->m_PlanBackwardMask );
  fftw_destroy_plan( this->m_PlanBackwardMask2 );

  fftw_destroy_plan( this->m_PlanFilterSquare );
  fftw_destroy_plan( this->m_PlanFilter );
  fftw_destroy_plan( this->m_PlanFilterMask );
  
  fftw_free( this->m_FilterMaskFT2 );
  fftw_free( this->m_FilterMaskFT );
  fftw_free( this->m_FilterSquareFT );
  fftw_free( this->m_FilterFT );
  fftw_free( this->m_ImageFT );
}

cmtk::TypedArray::SmartPtr
cmtk::SphereDetectionNormalizedBipolarMatchedFilterFFT::GetFilteredImageData( const Types::Coordinate sphereRadius, const int marginWidth )
{
  // check if the requested kernel parameters are the same we previously used.
  if ( (sphereRadius == this->m_SphereRadius) && (marginWidth == this->m_MarginWidth ) )
    return this->m_FilterResponse;

  this->m_SphereRadius = sphereRadius;
  this->m_MarginWidth = marginWidth;

  memset( this->m_FilterFT, 0, sizeof( fftw_complex ) * this->m_NumberOfPixels );
  memset( this->m_FilterSquareFT, 0, sizeof( fftw_complex ) * this->m_NumberOfPixels );
  memset( this->m_FilterMaskFT, 0, sizeof( fftw_complex ) * this->m_NumberOfPixels );

  this->MakeFilter( sphereRadius, marginWidth );

  const Types::DataItem denom2 = sqrt( this->m_SumFilterSquare - (this->m_SumFilter*this->m_SumFilter) / this->m_SumFilterMask );

  // compute filter kernel FT
  fftw_execute( this->m_PlanFilter );
  fftw_execute( this->m_PlanFilterSquare );
  fftw_execute( this->m_PlanFilterMask );
  
  // apply FT'ed filter to FT'ed image
  for ( size_t n = 0; n < this->m_NumberOfPixels; ++n )
    {
    this->m_FilterMaskFT2[n][0] = this->m_FilterMaskFT[n][0];
    this->m_FilterMaskFT2[n][1] = this->m_FilterMaskFT[n][1];

    FFTW::MultiplyInPlace( this->m_FilterMaskFT[n], this->m_ImageFT[n] );
    FFTW::MultiplyInPlace( this->m_FilterMaskFT2[n], this->m_ImageSquareFT[n] );
    FFTW::MultiplyInPlace( this->m_FilterFT[n], this->m_ImageFT[n] );
    }
  
  // transform filtered spectral data back into space domain
  fftw_execute( this->m_PlanBackward );
  fftw_execute( this->m_PlanBackwardMask );
  fftw_execute( this->m_PlanBackwardMask2 );

  for ( size_t n = 0; n < this->m_NumberOfPixels; ++n )
    {
    for ( int c = 0; c < 2; ++c )
      {
      this->m_FilterMaskFT[n][c] /= this->m_NumberOfPixels;
      this->m_FilterMaskFT2[n][c] /= this->m_NumberOfPixels;
      this->m_FilterFT[n][c] /= this->m_NumberOfPixels;
      }
    }
  
  // return filter response data
  if ( !this->m_FilterResponse )
    this->m_FilterResponse = TypedArray::Create( cmtk::TYPE_ITEM, this->m_NumberOfPixels );
  
  for ( size_t n = 0; n < this->m_NumberOfPixels; ++n )
    {
    const Types::DataItem num = FFTW::Magnitude( this->m_FilterFT[n] ) - (FFTW::Magnitude( this->m_FilterMaskFT[n] ) * this->m_SumFilter / this->m_SumFilterMask );
    const Types::DataItem denom1 = sqrt( std::max<Types::DataItem>( 0, FFTW::Magnitude( this->m_FilterMaskFT2[n] ) - FFTW::SumOfSquares( this->m_FilterMaskFT[n] ) / this->m_SumFilterMask ) );
    
    if ( (num == 0) || (denom1 ==0) )
      this->m_FilterResponse->Set( 0, n );
    else
      this->m_FilterResponse->Set( num / (denom1*denom2), n );
    }

  return this->m_FilterResponse;
}

void
cmtk::SphereDetectionNormalizedBipolarMatchedFilterFFT::MakeFilter( const Types::Coordinate sphereRadius, const int marginWidth )
{
  const int nRadius[3] = { 1 + marginWidth + static_cast<int>( sphereRadius / this->m_PixelSize[0] ), 
			   1 + marginWidth + static_cast<int>( sphereRadius / this->m_PixelSize[1] ), 
			   1 + marginWidth + static_cast<int>( sphereRadius / this->m_PixelSize[2] ) };

  this->m_SumFilter = this->m_SumFilterMask = this->m_SumFilterSquare = 0;
  
  // create a filter kernel for the sphere
  for ( int k = 0; k < nRadius[2]; ++k )
    for ( int j = 0; j < nRadius[1]; ++j )
      for ( int i = 0; i < nRadius[0]; ++i )
	{
	const cmtk::Types::Coordinate distance = sqrt( cmtk::MathUtil::Square( k * this->m_PixelSize[2] ) + cmtk::MathUtil::Square( j * this->m_PixelSize[1] ) + cmtk::MathUtil::Square( i * this->m_PixelSize[0] ) );
	if ( distance <= sphereRadius+marginWidth )
	  {
	  cmtk::Types::DataItem value = 1;
	  if ( (distance > sphereRadius) )
	    {
	    value = -1;
	    }

	  if ( value != 0 )
	    {
	    for ( int kk = k; kk < this->m_ImageDims[2]; kk += (this->m_ImageDims[2]-1-2*k) )
	      for ( int jj = j; jj < this->m_ImageDims[1]; jj += (this->m_ImageDims[1]-1-2*j) )
		for ( int ii = i; ii < this->m_ImageDims[0]; ii += (this->m_ImageDims[0]-1-2*i) )
		  {
		  const size_t ofs = ii+this->m_ImageDims[0] * (jj+this->m_ImageDims[1]*kk);
		  this->m_FilterFT[ofs][0] = value;
		  this->m_FilterSquareFT[ofs][0] = value*value;
		  this->m_FilterMaskFT[ofs][0] = 1;

		  this->m_SumFilter += value;
		  this->m_SumFilterSquare += value * value;
		  this->m_SumFilterMask += 1;
		  }
	    }
	  }
	}
}
