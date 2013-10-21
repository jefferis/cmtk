/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include "cmtkUniformVolumeMorphologicalOperators.h"

#include <Base/cmtkDistanceMap.h>
#include <Base/cmtkUniformDistanceMap.h>
#include <System/cmtkException.h>

#include <vector>

namespace
cmtk
{

UniformVolumeMorphologicalOperators
::UniformVolumeMorphologicalOperators( const UniformVolume::SmartConstPtr& uniformVolume ) 
  : m_UniformVolume( uniformVolume ) 
{}

TypedArray::SmartPtr
UniformVolumeMorphologicalOperators::GetErodedByDistance( const Types::Coordinate erodeBy ) const
{
  if ( !this->m_UniformVolume->GetData() )
    return TypedArray::SmartPtr( NULL );
  
  TypedArray::SmartPtr erodedData = UniformDistanceMap<Types::Coordinate>( *(this->m_UniformVolume), DistanceMap::INSIDE ).Get()->GetData();
  erodedData->Binarize( erodeBy + 0.5 );
  return erodedData->Convert( TYPE_BYTE );
}

TypedArray::SmartPtr
UniformVolumeMorphologicalOperators::GetDilatedByDistance( const Types::Coordinate dilateBy ) const
{
  if ( !this->m_UniformVolume->GetData() )
    return TypedArray::SmartPtr( NULL );
  
  TypedArray::SmartPtr dilatedData = UniformDistanceMap<Types::Coordinate>( *(this->m_UniformVolume) ).Get()->GetData();
  dilatedData->Binarize( dilateBy + 0.5 );
  dilatedData->Rescale( -1 /*scale*/, +1 /*offset*/ ); // this is binary inversion, 0->1, 1->0
  return dilatedData->Convert( TYPE_BYTE );
}


} // namespace cmtk
