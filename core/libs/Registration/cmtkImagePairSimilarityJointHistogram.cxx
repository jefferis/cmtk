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

#include "Registration/cmtkImagePairSimilarityJointHistogram.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairSimilarityJointHistogram::ImagePairSimilarityJointHistogram
( UniformVolume::SmartConstPtr& refVolume, UniformVolume::SmartConstPtr& fltVolume, const Interpolators::InterpolationEnum interpolation )
  : ImagePairSimilarityMeasure( Self::PrescaleData( refVolume, &this->m_NumberOfBinsX, &this->m_ScaleFactorReference, &this->m_ScaleOffsetReference ), 
				Self::PrescaleData( fltVolume, &this->m_NumberOfBinsY, &this->m_ScaleFactorFloating, &this->m_ScaleOffsetFloating ), interpolation )
{
  this->m_JointHistogram.Resize( this->m_NumberOfBinsX, this->m_NumberOfBinsY );
}

void
ImagePairSimilarityJointHistogram::SetReferenceVolume( const UniformVolume::SmartConstPtr& refVolume )
{
  Superclass::SetReferenceVolume( Self::PrescaleData( refVolume, &this->m_NumberOfBinsX, &this->m_ScaleFactorReference, &this->m_ScaleOffsetReference ) );
  this->m_JointHistogram.Resize( this->m_NumberOfBinsX, this->m_NumberOfBinsY );
}

void
ImagePairSimilarityJointHistogram::SetFloatingVolume( const UniformVolume::SmartConstPtr& fltVolume )
{
  Superclass::SetFloatingVolume( Self::PrescaleData( fltVolume, &this->m_NumberOfBinsY, &this->m_ScaleFactorFloating, &this->m_ScaleOffsetFloating ) );
  this->m_JointHistogram.Resize( this->m_NumberOfBinsX, this->m_NumberOfBinsY );
}

UniformVolume::SmartPtr
ImagePairSimilarityJointHistogram::PrescaleData
( const UniformVolume::SmartConstPtr& volume, size_t* numberOfBins, Types::DataItem* scaleFactor, Types::DataItem* scaleOffset )
{
  UniformVolume::SmartPtr newVolume( volume->CloneGrid() );
  newVolume->CreateDataArray( TYPE_ITEM );
  const size_t numberOfPixels = volume->GetNumberOfPixels();
    
  Types::DataItem value = 0;
  Types::DataItem minValue = FLT_MAX;
  Types::DataItem maxValue = -FLT_MAX;

  const DataGrid::IndexType& cropFrom = volume->CropRegion().From();
  const DataGrid::IndexType& cropTo = volume->CropRegion().To();
  const DataGrid::IndexType increments = volume->GetCropRegionIncrements();

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

    *scaleOffset = -minValue;
    *scaleFactor = 1.0;

    for ( size_t idx = 0; idx < numberOfPixels; ++idx ) 
      {
      if ( volume->GetDataAt( value, idx ) )
	newVolume->SetDataAt( static_cast<Types::DataItem>( value + *scaleOffset ), idx );
      else
	newVolume->GetData()->SetPaddingAt( idx );
      }
    }
    break;
    default: // Handle everything else as grey-level data.
    case DATACLASS_GREY: 
    {
    *numberOfBins = JointHistogramBase::CalcNumBins( volume );
    
    *scaleFactor = (*numberOfBins-1) / ( maxValue - minValue );
    *scaleOffset = -minValue * *scaleFactor;

    for ( size_t idx = 0; idx < numberOfPixels; ++idx ) 
      {
      if ( volume->GetDataAt( value, idx ) ) 
	{
	value = std::max( std::min( value, maxValue ), minValue );
	newVolume->SetDataAt( static_cast<Types::DataItem>( floor(*scaleFactor*value+*scaleOffset) ), idx );
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
