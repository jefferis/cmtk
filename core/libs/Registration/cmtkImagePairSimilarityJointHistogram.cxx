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

#include <cmtkImagePairSimilarityJointHistogram.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairSimilarityJointHistogram::ImagePairSimilarityJointHistogram
( const UniformVolume::SmartPtr& refVolume, const UniformVolume::SmartPtr& fltVolume, const Interpolators::InterpolationEnum interpolation )
  : ImagePairSimilarityMeasure( Self::PrescaleData( refVolume, &this->m_NumberOfBinsX ), Self::PrescaleData( fltVolume, &this->m_NumberOfBinsY ), interpolation )
{
  this->m_JointHistogram.Resize( this->m_NumberOfBinsX, this->m_NumberOfBinsY );
}

ImagePairSimilarityJointHistogram::ImagePairSimilarityJointHistogram
( const Self& other )
{
  StdErr << "Not implemented: " << __FILE__ << ":" << __LINE__ << "\n";
  exit(1);
}

ImagePairSimilarityJointHistogram::ImagePairSimilarityJointHistogram
( Self& other, const bool copyData )
  : ImagePairSimilarityMeasure( other, copyData )
{
  this->m_NumberOfBinsX = other.m_NumberOfBinsX;
  this->m_NumberOfBinsY = other.m_NumberOfBinsY;

  this->m_JointHistogram = other.m_JointHistogram;
}

UniformVolume::SmartPtr
ImagePairSimilarityJointHistogram::PrescaleData
( const UniformVolume::SmartPtr& volume, size_t* numberOfBins )
{
  UniformVolume::SmartPtr newVolume( volume->CloneGrid() );
  newVolume->CreateDataArray( TYPE_ITEM );
  const size_t numberOfPixels = volume->GetNumberOfPixels();
    
  Types::DataItem value = 0;
  Types::DataItem minValue = FLT_MAX;
  Types::DataItem maxValue = -FLT_MAX;

  int cropFrom[3], cropTo[3], increments[3];
  volume->GetCropRegion( cropFrom, cropTo, increments );
  int offset = increments[0];
  for ( int z = cropFrom[2]; z < cropTo[2]; ++z, offset += increments[2] ) 
    {
    for ( int y = cropFrom[1]; y < cropTo[1]; ++y, offset += increments[1] ) 
      {
      for ( int x = cropFrom[0]; x < cropTo[0]; ++x, ++offset ) 
	{
	if ( volume->GetDataAt( value, offset ) ) 
	  {
	  if ( value > maxValue ) maxValue = value;
	  if ( value < minValue ) minValue = value;
	  }
	}
      }
    }
  
  switch ( volume->GetData()->GetDataClass() ) 
    {
    case DATACLASS_LABEL: 
    {
    *numberOfBins = 1 + static_cast<unsigned int>(maxValue-minValue);
    if ( *numberOfBins > 254 ) 
      {
      StdErr << "Fatal error: Cannot handle more than 254 different labels.\n";
      exit( 1 );
      }
    
    for ( size_t idx = 0; idx < numberOfPixels; ++idx ) 
      {
      if ( volume->GetDataAt( value, idx ) )
	newVolume->SetDataAt( static_cast<Types::DataItem>( value - minValue ), idx );
      else
	newVolume->GetData()->SetPaddingAt( idx );
      }
    }
    break;
    default: // Handle everything else as grey-level data.
    case DATACLASS_GREY: 
    {
    *numberOfBins = JointHistogramBase::CalcNumBins( volume );
    const Types::DataItem binOffset = minValue;
    const Types::DataItem binWidth = ( maxValue - minValue ) / (*numberOfBins-1);
    const Types::DataItem factor = 1.0 / binWidth;
    
    for ( size_t idx = 0; idx < numberOfPixels; ++idx ) 
      {
      if ( volume->GetDataAt( value, idx ) ) 
	{
	value = std::max( std::min( value, maxValue ), minValue );
	newVolume->SetDataAt( static_cast<Types::DataItem>( floor(factor*(value-binOffset)) ), idx );
	} 
      else 
	{
	newVolume->GetData()->SetPaddingAt( idx );
	}
      }
    }
    break;
    } 

  return newVolume;
}

} // namespace cmtk
