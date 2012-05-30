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

#include "cmtkFitToXformListBase.h"

#include <Base/cmtkRegionIndexIterator.h>

cmtk::FitToXformListBase::FitToXformListBase( const UniformVolume& sampleGrid, const XformList& xformList, const bool absolute )
  : m_XformField( sampleGrid )
{
  this->m_XformValidAt.resize( this->m_XformField.GetNumberOfPixels() );
  std::fill( this->m_XformValidAt.begin(), this->m_XformValidAt.end(), true );

  size_t ofs = 0;
  for ( RegionIndexIterator<DataGrid::RegionType> voxelIt( this->m_XformField.GetWholeImageRegion() ); voxelIt != voxelIt.end(); ++voxelIt, ++ofs )
    {
    const Xform::SpaceVectorType v = this->m_XformField.GetGridLocation( voxelIt.Index() );
    
    Xform::SpaceVectorType u = v;
    if ( xformList.ApplyInPlace( u ) )
      {
      if ( !absolute )
	u -= v;
      this->m_XformField[ofs] = u;
      }
    else
      {
      this->m_XformValidAt[ofs] = false;
      }
    }
}
