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

#include "cmtkSphereDetectionBipolarMatchedFilterFFT.h"

cmtk::SphereDetectionBipolarMatchedFilterFFT::SphereDetectionBipolarMatchedFilterFFT( const UniformVolume& image )
  :   m_NumberOfPixels( image.GetNumberOfPixels() ),
      m_ImageDims( image.m_Dims ),
      m_PixelSize( image.m_Delta )
{
  this->m_ImageFT = fftw_alloc_complex( this->m_NumberOfPixels );
  this->m_FilterFT = fftw_alloc_complex( this->m_NumberOfPixels );

  this->m_PlanFilter = fftw_plan_dft_3d( this->m_ImageDims[2], this->m_ImageDims[1], this->m_ImageDims[0], this->m_FilterFT, this->m_FilterFT, FFTW_FORWARD, FFTW_ESTIMATE );
  this->m_PlanBackward = fftw_plan_dft_3d( this->m_ImageDims[2], this->m_ImageDims[1], this->m_ImageDims[0], this->m_FilterFT, this->m_FilterFT, FFTW_BACKWARD, FFTW_ESTIMATE );
  
  // initialize image FT
  fftw_plan plan_image = fftw_plan_dft_3d( this->m_ImageDims[2], this->m_ImageDims[1], this->m_ImageDims[0], this->m_ImageFT, this->m_ImageFT, FFTW_FORWARD, FFTW_ESTIMATE );
  for ( size_t n = 0; n < this->m_NumberOfPixels; ++n )
    {
    this->m_ImageFT[n][0] = image.GetDataAt( n );
    this->m_ImageFT[n][1] = 0;
    }
  
  fftw_execute( plan_image );
  fftw_destroy_plan( plan_image );
}

cmtk::SphereDetectionBipolarMatchedFilterFFT::~SphereDetectionBipolarMatchedFilterFFT()
{
  fftw_destroy_plan( this->m_PlanBackward );
  fftw_destroy_plan( this->m_PlanFilter );
  
  fftw_free( this->m_FilterFT );
  fftw_free( this->m_ImageFT );
}

inline void multiply( fftw_complex& lhs, const fftw_complex& rhs )
{
  fftw_complex product;
  product[0] = lhs[0]*rhs[0] - lhs[1]*rhs[1];
  product[1] = lhs[1]*rhs[0] + lhs[0]*rhs[1];

  lhs[0] = product[0];
  lhs[1] = product[1];
}

cmtk::TypedArray::SmartPtr
cmtk::SphereDetectionBipolarMatchedFilterFFT::GetFilteredImageData( const Types::Coordinate sphereRadius, const int marginWidth )
{
  memset( this->m_FilterFT, 0, sizeof( fftw_complex ) * this->m_NumberOfPixels );

  const size_t cnt = this->MakeFilter( sphereRadius, marginWidth );

  // compute filter kernel FT
  fftw_execute( this->m_PlanFilter );

  // apply FT'ed filter to FT'ed image
  for ( size_t n = 0; n < this->m_NumberOfPixels; ++n )
    {
    this->m_FilterFT[n][1] *= -1;  // correlation is multiplication with complex conjugate of filter, but we can also just conjugate the image
    multiply( this->m_FilterFT[n], this->m_ImageFT[n] );

    this->m_FilterFT[n][0] /= cnt;
    this->m_FilterFT[n][1] /= cnt;
    }
  
  // transform filtered spectral data back into space domain
  fftw_execute( this->m_PlanBackward );
  
  // return filtered magnitude image
  TypedArray::SmartPtr result = TypedArray::Create( cmtk::TYPE_ITEM, this->m_NumberOfPixels );
  for ( size_t n = 0; n < this->m_NumberOfPixels; ++n )
    {
    result->Set( sqrt( this->m_FilterFT[n][0]*this->m_FilterFT[n][0] + this->m_FilterFT[n][1]*this->m_FilterFT[n][1] ) / this->m_NumberOfPixels, n );
    }

  return result;
}

size_t
cmtk::SphereDetectionBipolarMatchedFilterFFT::MakeFilter( const Types::Coordinate sphereRadius, const int marginWidth )
{
  const int nRadius[3] = { 1 + marginWidth + static_cast<int>( sphereRadius / this->m_PixelSize[0] ), 
			   1 + marginWidth + static_cast<int>( sphereRadius / this->m_PixelSize[1] ), 
			   1 + marginWidth + static_cast<int>( sphereRadius / this->m_PixelSize[2] ) };

  // count number of non-zero filter elements for rough normalization
  size_t cnt = 0;

  // create a filter kernel for the sphere
  for ( int k = 0; k < nRadius[2]; ++k )
    for ( int j = 0; j < nRadius[1]; ++j )
      for ( int i = 0; i < nRadius[0]; ++i )
	{
	const cmtk::Types::Coordinate distance = sqrt( cmtk::MathUtil::Square( k * this->m_PixelSize[2] ) + cmtk::MathUtil::Square( j * this->m_PixelSize[1] ) + cmtk::MathUtil::Square( i * this->m_PixelSize[0] ) );
	if ( distance <= sphereRadius+marginWidth )
	  {
	  cmtk::Types::DataItem value = 0;
	  if ( distance >= sphereRadius-marginWidth )
	    {
	    value = 1;
	    }
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
		  this->m_FilterFT[ii+this->m_ImageDims[0] * (jj+this->m_ImageDims[1]*kk)][0] = value;
		  ++cnt;
		  }
	    }
	  }
	}
  
  return cnt;
}
