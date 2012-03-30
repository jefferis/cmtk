/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkUniformVolumeHoughTransformSphere.h"

#include <Base/cmtkRegionSphereSurfaceIterator.h>
#include <Base/cmtkRegionIndexIterator.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

TypedArray::SmartPtr
UniformVolumeHoughTransformSphere::Get( const Types::Coordinate radius ) const
{
  const UniformVolume& volume = *(this->m_UniformVolume);
  TypedArray::SmartPtr result( TypedArray::Create( volume.GetData()->GetType(), volume.GetNumberOfPixels() ) );

  const int dRadius[3] = { MathUtil::Round( radius / volume.m_Delta[0] ), MathUtil::Round( radius / volume.m_Delta[1] ), MathUtil::Round( radius / volume.m_Delta[2] ) };
  RegionSphereIterator<DataGrid::RegionType> sphereIterator( (DataGrid::IndexType( dRadius )) );

  const DataGrid::RegionType wholeImageRegion = volume.GetWholeImageRegion();

  const DataGrid::IndexType center = 0.5 * (wholeImageRegion.To() - wholeImageRegion.From());
  for ( sphereIterator = sphereIterator.begin(); sphereIterator != sphereIterator.end(); ++sphereIterator )
    {
    const DataGrid::IndexType pt = center+sphereIterator.Index();
    result->Set( 1, volume.GetOffsetFromIndex( pt ) );
    }

  return result;

#ifndef _OPENMP
  const DataGrid::RegionType region = wholeImageRegion;
#else // _OPENMP
  const int sliceFrom = wholeImageRegion.From()[2];
  const int sliceTo = wholeImageRegion.To()[2];
#pragma omp parallel for
  for ( int slice = sliceFrom; slice < sliceTo; ++slice )
    {
    DataGrid::RegionType region = wholeImageRegion;
    region.From()[2] = slice;
    region.To()[2] = slice+1;
#endif // #ifdef _OPENMP
    for ( RegionIndexIterator<DataGrid::RegionType> it( region ); it != it.end(); ++it )
      {
      const DataGrid::IndexType center = it.Index();
      const size_t toOffset = volume.GetOffsetFromIndex( center );
      for ( sphereIterator = sphereIterator.begin(); sphereIterator != sphereIterator.end(); ++sphereIterator )
	{
	const DataGrid::IndexType pt = center+sphereIterator.Index();
	if ( region.IsInside( pt ) )
	  {
	  result->Set( result->ValueAt( toOffset ) + volume.GetDataAt( volume.GetOffsetFromIndex( pt ) ), toOffset );
	  }
	}
      }
#ifdef _OPENMP
    }
#endif // #ifdef _OPENMP
  
  return result;
}

} // namespace cmtk
