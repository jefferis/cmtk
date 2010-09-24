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

#include "cmtkColormap.h"

#include <math.h>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#include <Base/cmtkTypedArray.h>
#include <Base/cmtkSegmentationLabel.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

const char *Colormap::StandardColormaps[] = 
{
  "Gray", "Red", "Green", "Blue", "Rainbow", "Labels", NULL
};

Colormap::Colormap()
{
  LookupTable = NULL;
  LookupTableEntries = 0;
  TableEntries = 256;
  SetStandardColormap( PALETTE_GRAY );
  DataRange[0] = 0;
  DataRange[1] = 4095;
  HueRange[0] = 0;
  HueRange[1] = 4095;
  ValueRange[0] = 0;
  ValueRange[1] = 4095;
  SaturationRange[0] = 0;
  SaturationRange[1] = 4095;
  Gamma = 0;
  Reverse = false;

  CreateSystemLabelColorMap( this->LabelColorMap );
}

Colormap::~Colormap()
{
  delete[] LookupTable;
}

void 
Colormap::SetFromStudy( const Study* study )
{
  if ( ! study ) return;

  // if there is a user-defined map, use that.
  if ( study->GetHaveUserColorMap() ) 
    {
    LabelColorMap = study->GetUserLabelMap();
    }
  
  // copy all other settings anyway, just in case.
  this->SetStandardColormap( study->GetStandardColormap() );
  this->SetReverse( study->GetReverseColormap() );
  this->SetDataRange( study->GetBlack(), study->GetWhite() );
  this->SetGamma( study->GetGamma() );
}

void
Colormap::SetStandardColormap( const int index )
{
  HaveUserMap = false;
  SetGamma( 0 );
  switch ( index ) 
    {
    case PALETTE_GRAY : 
      SetHueRange( 0, 0 );
      SetSaturationRange( 0, 0 );
      SetValueRange( 0, 1 );
      break;
    case PALETTE_RED : 
      SetHueRange( 0, 0 );
      SetSaturationRange( 1, 1 );
      SetValueRange( 0, 1 );
      break;
    case PALETTE_GREEN : 
      SetHueRange( 0.33, 0.33 );
      SetSaturationRange( 1, 1 );
      SetValueRange( 0, 1 );
      break;
    case PALETTE_BLUE : 
      SetHueRange( 0.66, 0.66 );
      SetSaturationRange( 1, 1 );
      SetValueRange( 0, 1 );
      break;
    case PALETTE_RAINBOW : 
      SetHueRange( 0.66, 0 );
      SetSaturationRange( 1, 1 );
      SetValueRange( 1, 1 );
      break;
    case PALETTE_LABELS:
    default:
      HaveUserMap = true;
      // nothing to do for user map; hopefully, there is one ;)
      break;
    }
}

void
Colormap::Apply( void *const outPtr, const TypedArray* inPtr, const bool generateAlpha )
{
  if ( ( outPtr == NULL ) || (inPtr == NULL) ) return;

  if ( (LookupTable == NULL) || (!TableEntries) ) 
    {
    memset( outPtr, 0, 3 * inPtr->GetDataSize() );
    return;
    }
  
  int size = inPtr->GetDataSize();
  if ( generateAlpha ) 
    {
    switch ( inPtr->GetType() ) 
      {
      case TYPE_BYTE:
	ApplyPrimitive<unsigned char>( (RGBA*) outPtr, (unsigned char*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((unsigned char*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_CHAR:
	ApplyPrimitive<char>( (RGBA*) outPtr, (char*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((char*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_SHORT:
	ApplyPrimitive<short>( (RGBA*) outPtr, (short*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((short*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_USHORT:
	ApplyPrimitive<unsigned short>( (RGBA*) outPtr, (unsigned short*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((unsigned short*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_INT:
	ApplyPrimitive<int>( (RGBA*) outPtr, (int*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((int*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_UINT:
	ApplyPrimitive<unsigned int>( (RGBA*) outPtr, (unsigned int*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((unsigned int*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_FLOAT:
	ApplyPrimitive<float>( (RGBA*) outPtr, (float*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((float*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_DOUBLE:
	ApplyPrimitive<double>( (RGBA*) outPtr, (double*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((double*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_NONE:
	break;
      }
    }
  else
    {
    switch ( inPtr->GetType() ) 
      {
      case TYPE_BYTE:
	ApplyPrimitive<unsigned char>( (RGB*) outPtr, (unsigned char*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((unsigned char*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_CHAR:
	ApplyPrimitive<char>( (RGB*) outPtr, (char*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((char*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_SHORT:
	ApplyPrimitive<short>( (RGB*) outPtr, (short*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((short*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_USHORT: 
	ApplyPrimitive<unsigned short>( (RGB*) outPtr, (unsigned short*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((unsigned short*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_INT:
	ApplyPrimitive<int>( (RGB*) outPtr, (int*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((int*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_UINT:
	ApplyPrimitive<unsigned int>( (RGB*) outPtr, (unsigned int*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((unsigned int*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_FLOAT:
	ApplyPrimitive<float>( (RGB*) outPtr, (float*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((float*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_DOUBLE:
	ApplyPrimitive<double>( (RGB*) outPtr, (double*)inPtr->GetDataPtr(), size, inPtr->GetPaddingFlag(), *((double*)inPtr->GetPaddingPtr()) );
	break;
      case TYPE_NONE:
	break;
      }
    }
}

template<class T>
void Colormap::ApplyPrimitive
( RGB *const outPtr, const T* inPtr, const unsigned int count,
  const bool paddingFlag, const T paddingData ) const
{
  if ( Reverse ) 
    {
    for ( unsigned int index = 0; index < count; ++index ) 
      {
      T data = inPtr[index];
      if ( (paddingFlag && (data == paddingData)) || !finite( data ) ) 
	data = 0;
      
      if ( data <= DataRange[0] ) 
	outPtr[index] = LookupTable[LookupTableEntries-1];
      else if ( data >= DataRange[1] ) 
	outPtr[index] = LookupTable[0];
      else
	outPtr[index] = LookupTable[ LookupTableEntries - 1 - (int)( (data - DataRange[0]) * InvDataRangeWidth ) ];
      }
    } 
  else
    {
    for ( unsigned int index = 0; index < count; ++index ) 
      {
      T data = inPtr[index];
      if ( (paddingFlag && (data == paddingData)) || !finite( data ) ) 
	data = 0;

      if ( data <= DataRange[0] ) 
	outPtr[index] = LookupTable[0];
      else if ( data >= DataRange[1] ) 
	outPtr[index] = LookupTable[LookupTableEntries-1];
      else
	outPtr[index] = LookupTable[ (int)( (data - DataRange[0]) * InvDataRangeWidth ) ];
      }
    }
}

template<class T>
void Colormap::ApplyPrimitive
( RGBA *const outPtr, const T* inPtr, const unsigned int count, const bool paddingFlag, const T paddingData ) const
{
  if ( Reverse ) 
    {
    for ( unsigned int index = 0; index < count; ++index ) 
      {
      T data = inPtr[index];
      if ( (paddingFlag && (data == paddingData)) || !finite( data ) ) 
	data = 0;
      
      if ( data <= DataRange[0] ) 
	outPtr[index] = LookupTable[LookupTableEntries-1];
      else if ( inPtr[index] >= DataRange[1] ) 
	outPtr[index] = LookupTable[0];
      else
	outPtr[index] = LookupTable[ LookupTableEntries - 1 - (int)( (data - DataRange[0]) * InvDataRangeWidth ) ];
      outPtr[index].Alpha = 255;
      }
    } 
  else
    {
    for ( unsigned int index = 0; index < count; ++index ) 
      {
      T data = inPtr[index];
      if ( (paddingFlag && (data == paddingData)) || !finite( data ) ) 
	data = 0;

      if ( data <= DataRange[0] ) 
	outPtr[index] = LookupTable[0];
      else if ( data >= DataRange[1] ) 
	outPtr[index] = LookupTable[LookupTableEntries-1];
      else
	outPtr[index] = LookupTable[ (int)( (data - DataRange[0]) * InvDataRangeWidth ) ];
      
      outPtr[index].Alpha = 255;
      }
    }
}

void Colormap::Execute ()
{
  if ( HaveUserMap ) 
    {
    // if user map exists, set table entry count to number of table entries.
    SegmentationLabelMap::const_iterator it = LabelColorMap.begin();
    int rangeFrom = it->first, rangeTo = it->first;
    while ( it != LabelColorMap.end() ) 
      {
      rangeFrom = std::min( rangeFrom, it->first );
      rangeTo = std::max( rangeTo, it->first );
      ++it;
      }

    TableEntries = (rangeTo - rangeFrom + 1);
    DataRange[0] = rangeFrom;
    DataRange[1] = rangeTo;
    } 
  else
    {
    TableEntries = 256;
    }
  
  if ( (LookupTable != NULL) && ( LookupTableEntries != TableEntries ) ) 
    {
    delete[] LookupTable;
    LookupTable = NULL;
    }
  
  if ( LookupTable == NULL ) 
    {
    LookupTable = Memory::AllocateArray<RGB>( TableEntries );
    LookupTableEntries = TableEntries;
    }
  
  if ( DataRange[0] != DataRange[1] )
    InvDataRangeWidth = (1.0 * (TableEntries - 1)) / (DataRange[1] - DataRange[0]);
  else
    InvDataRangeWidth = 0;
  
  if ( HaveUserMap ) 
    {
    // if user map exists, build table.
    for ( int index = 0; index < LookupTableEntries; ++index ) 
      {
      SegmentationLabelMap::const_iterator it = LabelColorMap.find( index );
      if ( it != LabelColorMap.end() ) 
	{
	const byte* rgb = it->second.GetRGB();
	LookupTable[index].R = rgb[0];
	LookupTable[index].G = rgb[1];
	LookupTable[index].B = rgb[2];
	} 
      else
	{
	LookupTable[index].R = LookupTable[index].G = LookupTable[index].B = 0;
	}
      }
    } 
  else 
    {
    // if no user-defined map, create map from HSV ramps.
    Types::DataItem H = HueRange[0];
    const Types::DataItem Hstep = (HueRange[1] - HueRange[0]) / (LookupTableEntries - 1);
    
    Types::DataItem S = SaturationRange[0];
    const Types::DataItem Sstep = (SaturationRange[1] - SaturationRange[0]) / (LookupTableEntries - 1);
    
    Types::DataItem V = ValueRange[0];
    const Types::DataItem Vstep = (ValueRange[1] - ValueRange[0]) / (LookupTableEntries - 1);
    
    if ( Gamma > 0 ) 
      {
      for ( int index = 0; index < LookupTableEntries; ++index, H += Hstep, S += Sstep, V += Vstep ) 
	{
	if ( V > 0 ) 
	  {
	  Types::DataItem Vgamma = exp( log(V) * (1/Gamma) );
	  HSV2RGB( LookupTable[index], H, S, Vgamma );
	  } 
	else
	  {
	  HSV2RGB( LookupTable[index], H, S, V );
	  }
	}
      } 
    else
      {
      for ( int index = 0; index < LookupTableEntries; ++index, H += Hstep, S += Sstep, V += Vstep ) 
	{
	HSV2RGB( LookupTable[index], H, S, V );
	}
      }
    }
}

void Colormap::HSV2RGB( RGB& rgb, Types::DataItem H, Types::DataItem S, Types::DataItem V )
{
  const Types::DataItem max = 1.0;
  const Types::DataItem third = 1.0 / 3.0;

  Types::DataItem R, G, B;
  // compute rgb assuming S = 1.0;
  if (H >= 0.0 && H <= third) 
    { // red -> green
    G = H/third;
    R = 1.0 - G;
    B = 0.0;
    } 
  else
    if (H >= third && H <= 2.0*third) 
      { // green -> blue
      B = (H - third)/third;
      G = 1.0 - B;
      R = 0.0;
      } 
    else
      { // blue -> red
      R = (H - 2.0 * third)/third;
      B = 1.0 - R;
      G = 0.0;
      }
  
  // add Saturation to the equation.
  S = S / max;
  
  (R *= S) += (1.0 - S);
  (G *= S) += (1.0 - S);
  (B *= S) += (1.0 - S);
  
  // Use value to get actual RGB 
  // normalize RGB first then apply value
  V = 3 * V / (R + G + B);
  R *= V;
  G *= V;
  B *= V;
  
  if (R > max) R = max;
  if (G > max) G = max;
  if (B > max) B = max;

  rgb.R = static_cast<unsigned char>( floor(255 * R) );
  rgb.G = static_cast<unsigned char>( floor(255 * G) );
  rgb.B = static_cast<unsigned char>( floor(255 * B) );
}

} // namespace cmtk
