/*
//
//  Copyright 2009-2010 SRI International
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

#include "cmtkImageOperationThreshold.h"

cmtk::UniformVolume::SmartPtr
cmtk::ImageOperationThreshold::Apply( cmtk::UniformVolume::SmartPtr& volume )
{
  cmtk::TypedArray::SmartPtr volumeData = volume->GetData();

  if ( this->m_Binarize )
    {
    volumeData->Binarize( this->m_Threshold );
    }
  else
    {
    cmtk::Types::DataItemRange range = volumeData->GetRange();
    
    if ( this->m_Above )
      range.m_UpperBound = this->m_Threshold;
    else
      range.m_LowerBound = this->m_Threshold;
    
    if ( this->m_ToPadding )
      volumeData->ThresholdToPadding( range );
    else
      volumeData->Threshold( range );
    }
  
  return volume;
}
