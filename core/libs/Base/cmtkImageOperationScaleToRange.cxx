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

#include <cmtkImageOperationScaleToRange.h>

#include <cmtkCommandLine.h>

void
cmtk::ImageOperationScaleToRange::New( const char* range )
{
  double rangeFrom, rangeTo;
  if ( 2 == sscanf( range, "%lf:%lf", &rangeFrom, &rangeTo ) )
    {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationScaleToRange( Types::DataItemRange( rangeFrom, rangeTo ) ) ) );
    }
  else
    {
    throw CommandLine::Exception( "Range must be given as two floating point numbers separated by ':', e.g., '0.5:1.0'" );
    }
}

cmtk::UniformVolume::SmartPtr
cmtk::ImageOperationScaleToRange::Apply( cmtk::UniformVolume::SmartPtr& volume )
{
  cmtk::TypedArray::SmartPtr volumeData = volume->GetData();
  volumeData->RescaleToRange( this->m_ToRange );
  return volume;
}
