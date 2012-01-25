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

cmtk::SplineWarpXform::SmartPtr 
cmtk::FitSplineWarpToDeformationField::Fit( const Types::Coordinate finalSpacing, const Types::Coordinate initialSpacing )
{
  // compute the start spacing of multi-level approximation by doubling final spacing until user-defined initial spacing is exceeded.
  Types::Coordinate spacing = finalSpacing;
  while ( spacing < initialSpacing )
    {
    spacing *= 2;
    }

  SplineWarpXform* splineWarp = new SplineWarpXform( this->m_DeformationField->m_Domain, spacing );
  splineWarp->RegisterVolumePoints( this->m_DeformationField->m_Dims, this->m_DeformationField->m_Spacing, this->m_DeformationField->m_Offset );

  for ( ; spacing >= initialSpacing; spacing /= 2 )
    {
    
    if ( spacing > initialSpacing )
      {
      splineWarp->Refine();
      splineWarp->RegisterVolumePoints( this->m_DeformationField->m_Dims, this->m_DeformationField->m_Spacing, this->m_DeformationField->m_Offset );
      }
    }
  
  return cmtk::SplineWarpXform::SmartPtr( splineWarp );
}
