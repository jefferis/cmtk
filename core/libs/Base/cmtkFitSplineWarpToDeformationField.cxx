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

#include "cmtkFitSplineWarpToDeformationField.h"

#include <Base/cmtkRegionIndexIterator.h>

cmtk::FitSplineWarpToDeformationField::FitSplineWarpToDeformationField( DeformationField::SmartConstPtr dfield ) 
  : m_DeformationField( dfield ),
    m_DeformationFieldFOV( dfield->m_Offset, dfield->m_Domain )
{  
}

cmtk::DataGrid::RegionType
cmtk::FitSplineWarpToDeformationField::GetDeformationGridRange( const UniformVolume::CoordinateRegionType& region ) const
{
  return cmtk::DataGrid::RegionType();
}

cmtk::SplineWarpXform::SmartPtr 
cmtk::FitSplineWarpToDeformationField::Fit( const Types::Coordinate finalSpacing, const Types::Coordinate initialSpacing )
{
  // compute the start spacing of multi-level approximation by doubling final spacing until user-defined initial spacing is exceeded.
  Types::Coordinate spacing = finalSpacing;
  while ( spacing < initialSpacing )
    {
    spacing *= 2;
    }

  // initialize B-spline transformation
  SplineWarpXform* splineWarp = new SplineWarpXform( this->m_DeformationField->m_Domain, spacing );
  splineWarp->RegisterVolumePoints( this->m_DeformationField->m_Dims, this->m_DeformationField->m_Spacing, this->m_DeformationField->m_Offset );

  // loop until final control point spacing
  for ( ; spacing >= initialSpacing; spacing /= 2 )
    {
    // loop over all control points
    for ( size_t cp = 0; cp < splineWarp->m_NumberOfControlPoints; ++cp )
      {
      // volume of influence for the current control point
      const DataGrid::RegionType voi = this->GetDeformationGridRange( splineWarp->GetVolumeOfInfluence( 3 * cp, this->m_DeformationFieldFOV, 0 /*fastMode=off*/ ) );
      
      // iterate over all voxels influenced by current control point.
      for ( RegionIndexIterator<DataGrid::RegionType> it( voi ); it != it.end(); ++it )
	{
	}
      }
    
    // refine control point grid if necessary for the next iteration
    if ( spacing > initialSpacing )
      {
      splineWarp->Refine();
      splineWarp->RegisterVolumePoints( this->m_DeformationField->m_Dims, this->m_DeformationField->m_Spacing, this->m_DeformationField->m_Offset );
      }
    }
  
  return cmtk::SplineWarpXform::SmartPtr( splineWarp );
}
