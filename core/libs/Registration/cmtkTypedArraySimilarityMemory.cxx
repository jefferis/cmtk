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
#include <cmtkTypedArraySimilarityMemory.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{
void
TypedArraySimilarityMemory::GetRangeX
( const TypedArray* array, const size_t numBins, 
  Types::DataItem& min, Types::DataItem& max )
{
  if ( ! this->ValidX )
    this->NumberBinsX = numBins;

  if ( ! this->ValidX || this->RepeatCheck ) {
    array->GetRange( min, max );
    if ( ! this->ValidX ) {
      this->RangeMinX = min;
      this->RangeMaxX = max;
      this->ValidX = true;
    } else if ( (min < this->RangeMinX) || (max > this->RangeMaxX ) ) {
      Types::DataItem binDelta = 
	(this->RangeMaxX - this->RangeMinX) / (this->NumberBinsX - 1);
      
      if ( min < this->RangeMinX ) {
	size_t addBins = 
	  1 + static_cast<size_t>( (this->RangeMinX - min) / binDelta);
	this->RangeMinX -= ( binDelta * addBins );
	this->NumberBinsY += addBins;
      }

      if ( max > this->RangeMaxX ) {
	size_t addBins = 
	  1 + static_cast<size_t>( (max - this->RangeMaxX) / binDelta);
	this->RangeMaxX += ( binDelta * addBins );
	this->NumberBinsY += addBins;
      }
    }
  }
  
  min = this->RangeMinX;
  max = this->RangeMaxX;
}

void
TypedArraySimilarityMemory::GetRangeY
( const TypedArray* array, const size_t numBins, 
  Types::DataItem& min, Types::DataItem& max )
{
  if ( ! this->ValidY )
    this->NumberBinsY = numBins;

  if ( ! this->ValidY || this->RepeatCheck ) {
    array->GetRange( min, max );
    if ( ! this->ValidY ) {
      this->RangeMinY = min;
      this->RangeMaxY = max;
      this->ValidY = true;
    } else if ( (min < this->RangeMinY) || (max > this->RangeMaxY ) ) {
      Types::DataItem binDelta = 
	(this->RangeMaxY - this->RangeMinY) / (this->NumberBinsY - 1);
      
      if ( min < this->RangeMinY ) {
	size_t addBins = 
	  1 + static_cast<size_t>( (this->RangeMinY - min) / binDelta);
	this->RangeMinY -= ( binDelta * addBins );
	this->NumberBinsY += addBins;
      }
      
      if ( max > this->RangeMaxY ) {
	size_t addBins = 
	  1 + static_cast<size_t>( (max - this->RangeMaxY) / binDelta);
	this->RangeMaxY += ( binDelta * addBins );
	this->NumberBinsY += addBins;
      }
    }
  }
  
  min = this->RangeMinY;
  max = this->RangeMaxY;
}

JointHistogram<unsigned int>*
TypedArraySimilarityMemory::CreateHistogram
( const TypedArray* array0, const TypedArray* array1 )
{
  unsigned int dataSize = array0->GetDataSize();

  size_t numBins = std::max<unsigned>
    ( std::min<unsigned>
      ( static_cast<unsigned>( sqrt( (float)dataSize ) ),
	this->MaxNumBins ), this->MinNumBins );
  
  Types::DataItem minX, maxX, minY, maxY;
  this->GetRangeX( array0, numBins, minX, maxX );
  this->GetRangeY( array1, numBins, minY, maxY );
  
  JointHistogram<unsigned int>* histogram =
    new JointHistogram<unsigned int>( this->NumberBinsX, this->NumberBinsY );

  histogram->SetRangeX( minX, maxX );
  histogram->SetRangeY( minY, maxY );

  return histogram;
}

} // namespace
