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

#include <cmtkFusionHV.h>

#include <cmtkRGB.h>
#include <cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

FusionHV::FusionHV() 
{ 
  Contrast = 1;
  ImageIndexV = 0;
};

void FusionHV::Execute()
{
  if ( !Input[0] || !Input[1] ) return;

  ImageRGB *output = this->GetOutput();
  output->CopyStructure( Input[0] );
  output->SetAlphaChannel( IMAGE_RGB );
  output->GetDataPtr( true /* forceAlloc */ );

  RGBA rgb1, rgb2;
  RGB outrgb;
  const double invSqrt3 = 1.0/sqrt((double)3);

  unsigned int outputPixels = output->GetNumPixels();
  for ( unsigned int idx = 0; idx < outputPixels; ++idx ) 
    {
    Input[  ImageIndexV]->GetPixel( rgb1, idx );
    Input[1-ImageIndexV]->GetPixel( rgb2, idx );
    
    const double brightness = 
      invSqrt3 * 
      sqrt( (double) ( MathUtil::Square( static_cast<unsigned short>(rgb1.R) ) + 
		       MathUtil::Square( static_cast<unsigned short>(rgb1.G) ) +
		       MathUtil::Square( static_cast<unsigned short>(rgb1.B) ) ) );
    
    const byte maxnorm = std::max( rgb2.R, std::max( rgb2.G, rgb2.B) );
    const double factor = static_cast<float>(255-Contrast*(255-brightness))/maxnorm;
    
    outrgb.R = (byte) (rgb2.R * factor);
    outrgb.G = (byte) (rgb2.G * factor);
    outrgb.B = (byte) (rgb2.B * factor);
    
    output->SetPixel( idx, outrgb );
    }
}

} // namespace cmtk
