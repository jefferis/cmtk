/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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
#include <cmtkTypedArraySimilarityMemory.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{
const Types::DataItemRange
TypedArraySimilarityMemory::GetRangeX
( const TypedArray* array, const size_t numBins )
{
  if ( ! this->ValidX )
    this->NumberBinsX = numBins;

  if ( ! this->ValidX || this->RepeatCheck ) 
    {
    const Types::DataItemRange range = array->GetRange();
    if ( ! this->ValidX ) 
      {
      this->RangeX.m_LowerBound = range.m_LowerBound;
      this->RangeX.m_UpperBound = range.m_UpperBound;
      this->ValidX = true;
      } 
    else
      if ( (range.m_LowerBound < this->RangeX.m_LowerBound) || (range.m_UpperBound > this->RangeX.m_UpperBound ) ) 
	{
	Types::DataItem binDelta = (this->RangeX.m_UpperBound - this->RangeX.m_LowerBound) / (this->NumberBinsX - 1);
	
	if ( range.m_LowerBound < this->RangeX.m_LowerBound ) 
	  {
	  const size_t addBins = 1 + static_cast<size_t>( (this->RangeX.m_LowerBound - range.m_LowerBound) / binDelta);
	  this->RangeX.m_LowerBound -= ( binDelta * addBins );
	  this->NumberBinsY += addBins;
	  }
	
	if ( range.m_UpperBound > this->RangeX.m_UpperBound ) 
	  {
	  const size_t addBins = 1 + static_cast<size_t>( (range.m_UpperBound - this->RangeX.m_UpperBound) / binDelta);
	  this->RangeX.m_UpperBound += ( binDelta * addBins );
	  this->NumberBinsY += addBins;
	  }
	}
    }
  
  return this->RangeX;
}

const Types::DataItemRange
TypedArraySimilarityMemory::GetRangeY
( const TypedArray* array, const size_t numBins )
{
  if ( ! this->ValidY )
    this->NumberBinsY = numBins;

  if ( ! this->ValidY || this->RepeatCheck ) 
    {
    const Types::DataItemRange range = array->GetRange();
    if ( ! this->ValidY ) 
      {
      this->RangeY.m_LowerBound = range.m_LowerBound;
      this->RangeY.m_UpperBound = range.m_UpperBound;
      this->ValidY = true;
      } 
    else
      if ( (range.m_LowerBound < this->RangeY.m_LowerBound) || (range.m_UpperBound > this->RangeY.m_UpperBound ) ) 
	{
	Types::DataItem binDelta = (this->RangeY.m_UpperBound - this->RangeY.m_LowerBound) / (this->NumberBinsY - 1);
	
	if ( range.m_LowerBound < this->RangeY.m_LowerBound ) 
	  {
	  const size_t addBins = 1 + static_cast<size_t>( (this->RangeY.m_LowerBound - range.m_LowerBound) / binDelta);
	  this->RangeY.m_LowerBound -= ( binDelta * addBins );
	  this->NumberBinsY += addBins;
	  }
	
	if ( range.m_UpperBound > this->RangeY.m_UpperBound ) 
	  {
	  const size_t addBins = 1 + static_cast<size_t>( (range.m_UpperBound - this->RangeY.m_UpperBound) / binDelta);
	  this->RangeY.m_UpperBound += ( binDelta * addBins );
	  this->NumberBinsY += addBins;
	  }
	}
    }
  
  return this->RangeY;
}

JointHistogram<unsigned int>::SmartPtr
TypedArraySimilarityMemory::CreateHistogram
( const TypedArray* array0, const TypedArray* array1 )
{
  const unsigned int dataSize = array0->GetDataSize();
  const size_t numBins = std::max<unsigned>( std::min<unsigned>( static_cast<unsigned>( sqrt( (float)dataSize ) ), this->MaxNumBins ), this->MinNumBins );
  
  Types::DataItemRange rangeX = this->GetRangeX( array0, numBins );
  Types::DataItemRange rangeY = this->GetRangeY( array1, numBins );

  JointHistogram<unsigned int>::SmartPtr histogram( new JointHistogram<unsigned int>( this->NumberBinsX, this->NumberBinsY ) );
  
  histogram->SetRangeX( rangeX );
  histogram->SetRangeY( rangeY );

  return histogram;
}

} // namespace
