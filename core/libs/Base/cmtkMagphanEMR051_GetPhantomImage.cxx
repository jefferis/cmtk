/*
//
//  Copyright 2012 SRI International
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

#include "cmtkMagphanEMR051.h"

#include <Base/cmtkUniformVolumePainter.h>
#include <Base/cmtkAnatomicalOrientation.h>

cmtk::UniformVolume::SmartPtr
cmtk::MagphanEMR051::GetPhantomImage( const cmtk::Types::Coordinate resolution, const bool labels )
{
  const int npx = 1 + static_cast<int>( 200.0 / resolution );
  const int dims[3] = { npx, npx, npx };
  UniformVolume::SmartPtr result( new UniformVolume( DataGrid::IndexType::FromPointer( dims ), resolution, resolution, resolution ) );
  result->SetMetaInfo( cmtk::META_SPACE, cmtk::AnatomicalOrientation::ORIENTATION_STANDARD );
  result->SetMetaInfo( cmtk::META_SPACE_ORIGINAL, cmtk::AnatomicalOrientation::ORIENTATION_STANDARD );
  result->CreateDataArray( TYPE_SHORT );
  
  const Types::Coordinate offset[3] = { -100, -100, -100 };
  result->m_Offset = UniformVolume::CoordinateVectorType::FromPointer( offset );
  
  UniformVolumePainter painter( result, UniformVolumePainter::COORDINATES_ABSOLUTE );
  for ( int idx = 0; idx < 165; ++idx )
    {
    const Types::DataItem value = ( labels ) ? (idx+1) : Self::SphereTable[idx].m_EstimatedT1;
    painter.DrawSphere( UniformVolume::CoordinateVectorType::FromPointer( Self::SphereTable[idx].m_CenterLocation ), Self::SphereTable[idx].m_Diameter / 2, value );
    }
  
  return result;
}
