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

#include <cmtkFusionOverlay.h>

#include <cmtkRGB.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

FusionOverlay::FusionOverlay() 
{ 
  TopImageIndex = 0; 
  OriginalImage[0] = OriginalImage[1] = NULL;
};

FusionOverlay::~FusionOverlay()
{
  if ( OriginalImage[0] ) OriginalImage[0]->Delete();
  if ( OriginalImage[1] ) OriginalImage[1]->Delete();
}

void FusionOverlay::SetOriginalImage( const int index, Image *const originalImage )
{
  this->ReplaceObject( OriginalImage[index], originalImage );
}

void FusionOverlay::Execute()
{
  if ( !Input[0] || !Input[1] ) return;
  if ( !OriginalImage[0] || !OriginalImage[1] ) return;

  ImageRGB *output = this->GetOutput();
  output->CopyStructure( Input[0] );
  output->SetAlphaChannel( IMAGE_RGB );
  output->GetDataPtr( true /* forceAlloc */ );

  RGB rgb;
  unsigned int outputPixels = output->GetNumPixels();
  for ( unsigned int idx = 0; idx < outputPixels; ++idx ) 
    {
    const double topValue = OriginalImage[TopImageIndex]->GetDataAt( idx );
    
    if ( (topValue >= LowerThreshold) && (topValue <= UpperThreshold) )
      Input[  TopImageIndex]->GetPixel( rgb, idx );
    else
      Input[1-TopImageIndex]->GetPixel( rgb, idx );
    
    output->SetPixel( idx, rgb );
  }
}

} // namespace cmtk
