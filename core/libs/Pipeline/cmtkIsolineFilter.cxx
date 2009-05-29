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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkIsolineFilter.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

void IsolineFilter::Execute()
{
  if ( Input == NULL ) return;
  if ( ! NumberOfLevels ) return;
  
  const TypedArray *inData = Input->GetData();

  Image *output = this->GetOutput();
  output->CopyStructure( Input );
  output->SetDataType( Input->GetDataType() );
  TypedArray::SmartPtr outData = output->GetData();
  outData->ClearArray();

  unsigned int width, height;
  Input->GetDims( width, height );

  Types::DataItem item00, item01, item10, item11;

  float rangeIncrement = (RangeTo - RangeFrom );
  if ( NumberOfLevels > 1 )
    rangeIncrement /= (NumberOfLevels-1);

  unsigned int offset = 0;
  for ( unsigned int y=0; y<height-1; ++y, ++offset ) {
    inData->Get( item10, offset );
    inData->Get( item11, offset+width );
    for ( unsigned int x=0; x<width-1; ++x, ++offset ) {
      item00 = item10;
      item01 = item11;
      inData->Get( item10, offset+1 );
      inData->Get( item11, offset+width+1 );

      float levelValue = RangeFrom;
      for ( unsigned int level=0; level<NumberOfLevels; 
	    ++level, levelValue+=rangeIncrement ) {
	if ( (item00 >= levelValue) || (item01 >= levelValue) ||
	     (item10 >= levelValue) || (item11 >= levelValue) ) {
	  if ( item00 < levelValue ) 
	    outData->Set( levelValue, offset );
	  if ( item10 < levelValue ) 
	    outData->Set( levelValue, offset+1 );
	  if ( item01 < levelValue ) 
	    outData->Set( levelValue, offset+width );
	  if ( item11 < levelValue ) 
	    outData->Set( levelValue, offset+width+1 );
	}
      }
    }
  }
  
  this->UpdateExecuteTime();
}

} // namespace cmtk
