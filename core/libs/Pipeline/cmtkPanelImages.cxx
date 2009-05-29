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

#include <cmtkPanelImages.h>

#include <cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

void PanelImages::Execute()
{
  if ( (Input[0] == NULL) || (Input[1] == NULL) ) return;

  const byte *data0 = Input[0]->GetDataPtr();
  const byte *data1 = Input[1]->GetDataPtr();

  if ( (data0 == NULL) || (data1 == NULL) ) return;

  ImageRGB *output = this->GetOutput();

  unsigned int Dims0[2], Dims1[2];
  Input[0]->GetDims( Dims0 );
  Input[1]->GetDims( Dims1 );
  unsigned int maxY = std::max( Dims0[1], Dims1[1] );
  
  output->SetDims( Dims0[0] + Dims1[0], maxY );
  output->SetSpacing( Input[0]->GetSpacing() );
  output->SetAlphaChannel( IMAGE_RGB );
  
  byte *dataDst = output->GetDataPtr( true /*forceAlloc*/ );

  unsigned int offsetDst = 0, offset0 = 0, offset1 = 0;
  for ( unsigned int y = 0; y < maxY; ++y ) 
    {
    if ( y < Dims0[1] ) 
      {
      memcpy( dataDst + offsetDst, data0 + offset0, 3 * Dims0[0] );
      offset0 += 3 * Dims0[0];
      } 
    else
      {
      memset( dataDst + offsetDst, 0, 3 * Dims0[0] );
      }
    offsetDst += 3 * Dims0[0];
    if ( y < Dims1[1] ) 
      {
      memcpy( dataDst + offsetDst, data1 + offset1, 3 * Dims1[0] );
      offset1 += 3 * Dims1[0];
      } 
    else
      {
      memset( dataDst + offsetDst, 0, 3 * Dims1[0] );
      }
    offsetDst += 3 * Dims1[0];
    }
  
  this->UpdateExecuteTime();
}

} // namespace cmtk
