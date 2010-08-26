/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "Pipeline/cmtkImageToImageRGB.h"

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

ImageToImageRGB::ImageToImageRGB() :
  CheckerboxPadding( true )
{
  this->m_Image = NULL;
  this->m_Colormap = NULL;
  AlphaMode = AlphaModeNone;

  this->RegisterInput( &this->m_Image );
  this->RegisterInput( &this->m_Colormap );
}

ImageToImageRGB::~ImageToImageRGB()
{
  if ( this->m_Image ) 
    this->m_Image->Unregister();
  if ( this->m_Colormap ) 
    this->m_Colormap->Unregister();
}

void ImageToImageRGB::SetInput( Image *const image )
{
  this->ReplaceObject( this->m_Image, image );
}

void ImageToImageRGB::SetColormap( Colormap *const colormap )
{
  this->ReplaceObject( this->m_Colormap, colormap );
}

void ImageToImageRGB::Execute()
{
  if ( (this->m_Image == NULL) || (this->m_Colormap == NULL) ) return;
  const TypedArray *inPtr = this->m_Image->GetData(); 
  if ( !inPtr ) return;

  ImageRGB *output = this->GetOutput();
  output->CopyStructure( this->m_Image );

  if ( AlphaMode == AlphaModeNone )
    output->SetAlphaChannel( IMAGE_RGB );
  else
    output->SetAlphaChannel( IMAGE_RGBA );
  
  void *outPtr = output->GetDataPtr( true /* forceAlloc */ );

  switch ( AlphaMode ) 
    {
    case AlphaModeNone:
      this->m_Colormap->Apply( outPtr, inPtr );
      if ( inPtr->GetPaddingFlag() )
	this->MarkPaddingData( output->GetDims( AXIS_X ), output->GetDims( AXIS_Y ), static_cast<RGB*>( outPtr ), inPtr );
      break;
    case AlphaModeConst:
      this->m_Colormap->Apply( outPtr, inPtr, true /* generateAlpha */ );
      if ( inPtr->GetPaddingFlag() )
	this->MarkPaddingData( output->GetDims( AXIS_X ), output->GetDims( AXIS_Y ), static_cast<RGBA*>( outPtr ), inPtr );
      break;
    }
  
  this->UpdateExecuteTime();
}

template<class T> 
void
ImageToImageRGB::MarkPaddingData
( const unsigned int dimsx, const unsigned int dimsy, T *const rgba, const TypedArray* data ) const
{
  T* p = rgba;
  unsigned int idx = 0;

  byte bright = 170;
  byte dark = 80;
  if ( !this->CheckerboxPadding )
    bright = dark = 0;
  
  for ( unsigned int y = 0; y < dimsy; ++y ) 
    {
    for ( unsigned int x = 0; x < dimsx; ++x, ++idx, ++p ) 
      {
      if ( data->PaddingDataAt( idx ) ) 
	{
	if ( ((x >> 4)&1) ^ ((y >> 4)&1) ) 
	  {
	  p->R = p->G = p->B = bright;
	  } 
	else
	  {
	  p->R = p->G = p->B = dark;
	  }
	}
      }
    }
}

} // namespace cmtk
