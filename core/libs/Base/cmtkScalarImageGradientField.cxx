/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include "cmtkScalarImageGradientField.h"

#include <Base/cmtkRegionIndexIterator.h>

cmtk::ScalarImageGradientField::ScalarImageGradientField( const UniformVolume& volume )
  : m_GradientField( new Self::GradientFieldType( volume.m_Dims, volume.m_Size ) )
{
  const DataGrid::RegionType wholeImageRegion = volume.GetWholeImageRegion();

  size_t ofs = 0;
  for ( RegionIndexIterator<DataGrid::RegionType> it( wholeImageRegion ); it != it.end(); ++it, ++ofs )
    {
    const DataGrid::IndexType idx = it.Index();
    
    for ( int dim = 0; dim < 3; ++dim )
      {
      Types::Coordinate div = 0;
      
      const DataGrid::IndexType idxUp = idx.AddScalarToOne( dim, 1 );
      if ( idxUp[dim] < wholeImageRegion.To()[dim] )
	{
	(*this->m_GradientField)[ofs][dim] = volume.GetDataAt( volume.GetOffsetFromIndex( idxUp ) );
	div += 1.0;
	}
      else
	(*this->m_GradientField)[ofs][dim] = volume.GetDataAt( volume.GetOffsetFromIndex( idx ) );
      
      const DataGrid::IndexType idxDown = idx.AddScalarToOne( dim, -1 );
      if ( idx[dim] > wholeImageRegion.From()[dim] )
	{
	(*this->m_GradientField)[ofs][dim] -= volume.GetDataAt( volume.GetOffsetFromIndex( idxDown ) );
	div += 1.0;
	}
      else
	(*this->m_GradientField)[ofs][dim] -= volume.GetDataAt( volume.GetOffsetFromIndex( idx ) );
      
      (*this->m_GradientField)[ofs][dim] /= div;
      }
    }
}

