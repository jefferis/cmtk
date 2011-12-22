/*
//
//  Copyright 2009-2011 SRI International
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

#include "cmtkImageOperationRegionFilter.h"

void
cmtk::ImageOperationRegionFilter
::NewGeneric( const Self::OperatorEnum op, const char* arg )
{
  int radiusX = 1;
  int radiusY = 1;
  int radiusZ = 1;
  
  const size_t nRadii = sscanf( arg, "%10d,%10d,%10d", &radiusX, &radiusY, &radiusZ );
  if ( nRadii == 1 )
    {
    radiusZ = radiusY = radiusX;
    }
  else
    {
    if ( nRadii != 3 )
      {
      cmtk::StdErr << "ERROR: downsampling radii must either be three integers, x,y,z, or a single integer\n";
      exit( 1 );
      }
    }
  ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationRegionFilter( op, radiusX, radiusY, radiusZ ) ) );
}

cmtk::UniformVolume::SmartPtr
cmtk::ImageOperationRegionFilter
::Apply( cmtk::UniformVolume::SmartPtr& volume )
{
  switch ( this->m_Operator )
    {
    case MEDIAN:
      volume->SetData( DataGridFilter( volume ).GetDataMedianFiltered( this->m_RadiusX, this->m_RadiusY, this->m_RadiusZ ) );
      break;
    case MEAN:
      volume->SetData( DataGridFilter( volume ).RegionMeanFilter( this->m_RadiusX, this->m_RadiusY, this->m_RadiusZ ) );
      break;
    case FAST_MEAN:
      volume->SetData( DataGridFilter( volume ).FastRegionMeanFilter( this->m_RadiusX, this->m_RadiusY, this->m_RadiusZ ) );
      break;
    case VARIANCE:
      volume->SetData( DataGridFilter( volume ).RegionVarianceFilter( this->m_RadiusX, this->m_RadiusY, this->m_RadiusZ ) );
      break;
    case FAST_VARIANCE:
      volume->SetData( DataGridFilter( volume ).FastRegionVarianceFilter( this->m_RadiusX, this->m_RadiusY, this->m_RadiusZ ) );
      break;
    case THIRD_MOMENT:
      volume->SetData( DataGridFilter( volume ).RegionThirdMomentFilter( this->m_RadiusX, this->m_RadiusY, this->m_RadiusZ ) ); 
      break;
    case STANDARD_DEVIATION:
      volume->SetData( DataGridFilter( volume ).RegionStandardDeviationFilter( this->m_RadiusX, this->m_RadiusY, this->m_RadiusZ ) );
      break;
    case SMOOTHNESS:
      volume->SetData( DataGridFilter( volume ).RegionSmoothnessFilter( this->m_RadiusX, this->m_RadiusY, this->m_RadiusZ ) );
      break;
    }
  return volume;
}
