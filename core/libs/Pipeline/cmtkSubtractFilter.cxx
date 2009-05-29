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

#include <cmtkSubtractFilter.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

void
SubtractFilter::Execute()
{
  if ( (Input[0] == NULL) || (Input[1] == NULL) ) return;

  const TypedArray* inData0 = ( SubtractImageIndex == 0 ) ? Input[0]->GetData() : Input[1]->GetData();
  const TypedArray* inData1 = ( SubtractImageIndex == 0 ) ? Input[1]->GetData() : Input[0]->GetData();

  Image *output = this->GetOutput();
  output->CopyStructure( Input[0] );
  //  output->SetDataType( Input[0]->GetDataType() );
  output->SetDataType( TYPE_FLOAT );
  TypedArray::SmartPtr outData = output->GetData();

  unsigned int width, height;
  Input[0]->GetDims( width, height );
  unsigned int size = width * height;

  Types::DataItem item0, item1;
  if ( Absolute ) 
    {
    for ( unsigned int offset=0; offset<size; ++offset ) 
      {
      inData0->Get( item0, offset );
      inData1->Get( item1, offset );
      if ( UseThresholds && ( ( item0 < Thresholds[0] ) || ( item0 > Thresholds[1] ) || ( item1 < Thresholds[0] ) || ( item1 > Thresholds[1] ) ) )
	outData->Set( 0, offset );
      else
	outData->Set( fabs( item0 - item1 ), offset );
      }
    } 
  else
    {
    for ( unsigned int offset=0; offset<size; ++offset ) 
      {
      inData0->Get( item0, offset );
      inData1->Get( item1, offset );
      if ( UseThresholds && ( ( item0 < Thresholds[0] ) || ( item0 > Thresholds[1] ) || ( item1 < Thresholds[0] ) || ( item1 > Thresholds[1] ) ) )
	outData->Set( 0, offset );
      else
	outData->Set( item0 - item1, offset );
      }
    }
  
  this->UpdateExecuteTime();
}

} // namespace cmtk
