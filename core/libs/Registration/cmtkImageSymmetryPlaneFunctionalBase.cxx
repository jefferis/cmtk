/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include "cmtkImageSymmetryPlaneFunctionalBase.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImageSymmetryPlaneFunctionalBase::ImageSymmetryPlaneFunctionalBase
( UniformVolume::SmartConstPtr& volume ) 
  : m_Volume( volume ),
    m_FixOffset( false )
{
}

ImageSymmetryPlaneFunctionalBase::ImageSymmetryPlaneFunctionalBase
( UniformVolume::SmartConstPtr& volume, 
  const Types::DataItemRange& valueRange )
  : m_Volume( Self::ApplyThresholds( *volume, valueRange ) ),
    m_FixOffset( false )
{
}

Types::Coordinate 
ImageSymmetryPlaneFunctionalBase::GetParamStep 
( const size_t idx, const Types::Coordinate mmStep ) 
  const
{
  switch ( idx ) 
    {
    // plane offset is a translation
    case 0:
      if ( this->m_FixOffset )
	return 0;
      else
	return mmStep;
      // the other two parameters are rotations
    case 1:
    case 2:
      return mmStep / sqrt( MathUtil::Square( 0.5 * m_Volume->m_Size[0] ) + MathUtil::Square( 0.5 * m_Volume->m_Size[1] ) + MathUtil::Square( 0.5 * m_Volume->m_Size[2] ) ) * 90/M_PI;
    }
  return mmStep;
}

UniformVolume::SmartPtr 
ImageSymmetryPlaneFunctionalBase::ApplyThresholds( const UniformVolume& volume, const Types::DataItemRange& valueRange )
{
  UniformVolume::SmartPtr result( volume.Clone() );
  result->GetData()->Threshold( valueRange );
  return result;
}


} // namespace cmtk
