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

  size_t ofsPlusMinus = 1;
  for ( int dim = 0; dim < 3; ++dim )
    {
    size_t ofs = 0;
    for ( RegionIndexIterator<DataGrid::RegionType> it( wholeImageRegion ); it != it.end(); ++it, ++ofs )
      {
      const DataGrid::IndexType idx = it.Index();
      Types::Coordinate div = 0;
      
      if ( idx[dim]+1 < wholeImageRegion.To()[dim] )
	{
	(*this->m_GradientField)[ofs][dim] = volume.GetDataAt( ofs + ofsPlusMinus );
	div += 1.0;
	}
      else
	(*this->m_GradientField)[ofs][dim] = volume.GetDataAt( ofs );
      
      if ( idx[dim]-1 > wholeImageRegion.From()[dim] )
	{
	(*this->m_GradientField)[ofs][dim] -= volume.GetDataAt( ofs - ofsPlusMinus );
	div += 1.0;
	}
      else
	(*this->m_GradientField)[ofs][dim] -= volume.GetDataAt( ofs );
      
      (*this->m_GradientField)[ofs][dim] /= div;
      }
    
    // compute increment for next dimension
    ofsPlusMinus *= volume.m_Dims[dim];
    }
}

