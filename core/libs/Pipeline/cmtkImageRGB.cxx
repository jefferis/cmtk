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

#include <cmtkconfig.h>

#include <cmtkImageRGB.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

ImageRGB *ImageRGB::New()
{
  return new ImageRGB;
}

ImageRGB::ImageRGB()
{
  Data = NULL;
  DataSize = 0;
  AlphaChannel = IMAGE_RGB;
  BytesPerPixel = 3;
}

ImageRGB::~ImageRGB()
{
  if ( Data ) delete[] Data;
}

byte* 
ImageRGB::GetDataPtr( const bool forceAlloc )
{
  if ( ! forceAlloc ) return Data;

  if ( Data == NULL ) 
    {
    DataSize = BytesPerPixel * this->Plane::GetNumPixels();
    Data = Memory::AllocateArray<byte>( DataSize );
    } 
  else
    {
    if ( DataSize != (BytesPerPixel * this->Plane::GetNumPixels()) ) 
      {
      delete[] Data;
      Data = NULL;
      return this->GetDataPtr( true /* forceAlloc */ );
      } 
    }
  return Data;
}

void ImageRGB::SetAlphaChannel( const ImageAlphaToggle alphaChannel, const bool convertData )
{
  if ( alphaChannel != AlphaChannel ) 
    {
    AlphaChannel = alphaChannel;
    BytesPerPixel = (AlphaChannel == IMAGE_RGB) ? 3 : 4;
    
    byte *oldData = Data;
    // fake initialization by setting Data to NULL
    Data = NULL;
    this->GetDataPtr( true /* forceAlloc */ );
    
    // convert old image data
    if ( convertData ) 
      {
      byte *fromPtr = oldData;
      byte *toPtr = Data;
      unsigned int numberOfPixels = this->GetNumPixels();
      
      if ( AlphaChannel == IMAGE_RGB ) 
	{
	// from RGBA to RGB: remove alpha value
	for ( unsigned int i = 0; i < numberOfPixels; ++i, fromPtr += 4, toPtr += 3 ) 
	  {
	  toPtr[0] = fromPtr[0];
	  toPtr[1] = fromPtr[1];
	  toPtr[2] = fromPtr[2];
	  }
	} 
      else
	{
	// from RGB to RGBA: set alpha to opaque
	for ( unsigned int i = 0; i < numberOfPixels; ++i, fromPtr += 3, toPtr += 4 ) 
	  {
	  toPtr[0] = fromPtr[0];
	  toPtr[1] = fromPtr[1];
	  toPtr[2] = fromPtr[2];
	  toPtr[3] = 255;
	  }
	}
      }
    
    delete[] oldData;
    }
}

void ImageRGB::GetPixel( RGBA& rgb, const int index )
{
  byte* pixelPtr = Data + (BytesPerPixel * index);
  rgb.R = pixelPtr[0];
  rgb.G = pixelPtr[1];
  rgb.B = pixelPtr[2];
  if ( AlphaChannel == IMAGE_RGBA )
    rgb.Alpha = pixelPtr[3];
  else
    rgb.Alpha = 255;
}

void ImageRGB::SetPixel( const int index, const RGBA& rgb )
{
  byte* pixelPtr = Data + (BytesPerPixel * index);
  pixelPtr[0] = rgb.R;
  pixelPtr[1] = rgb.G;
  pixelPtr[2] = rgb.B;
  if ( AlphaChannel == IMAGE_RGBA )
    pixelPtr[3] = rgb.Alpha;
}

void ImageRGB::GetPixel( RGB& rgb, const int index )
{
  byte* pixelPtr = Data + (BytesPerPixel * index);
  rgb.R = pixelPtr[0];
  rgb.G = pixelPtr[1];
  rgb.B = pixelPtr[2];
}

void ImageRGB::SetPixel( const int index, const RGB& rgb )
{
  byte* pixelPtr = Data + (BytesPerPixel * index);
  pixelPtr[0] = rgb.R;
  pixelPtr[1] = rgb.G;
  pixelPtr[2] = rgb.B;
}

bool 
ImageRGB::IsGreyscale() const
{
  unsigned int numPixels = this->GetNumPixels();
  const byte *data = this->GetDataPtr();
  byte increment = (this->GetAlphaChannel() == IMAGE_RGB) ? 3 : 4;

  for ( unsigned int idx = 0; idx < numPixels; ++idx, data += increment ) 
    {
    if ( (data[0] != data[1]) || (data[1] != data[2]) )
      return 0;
    }
  
  return 1;
}

} // namespace cmtk
