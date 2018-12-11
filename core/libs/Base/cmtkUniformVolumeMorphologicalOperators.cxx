/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2014 SRI International
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
#include <set>

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
  erodedData->SetDataClass( DATACLASS_LABEL );

  return erodedData->Convert( TYPE_BYTE );
}

TypedArray::SmartPtr
UniformVolumeMorphologicalOperators::GetErodedByDistanceMultiLabels( const Types::Coordinate erodeBy ) const
{
  if ( !this->m_UniformVolume->GetData() )
    return TypedArray::SmartPtr( NULL );

  const UniformVolume& volume = *(this->m_UniformVolume);
  const size_t nPixels = volume.GetNumberOfPixels();

  // Make a set of label values that occur in the volume. Treat everything as unsigned int at most.
  unsigned int maxLabel = 0;
  std::set<unsigned int> existingLabels;
  for ( size_t idx = 0; idx < nPixels; ++idx )
    {
    const unsigned int value = static_cast<int>( volume.GetDataAt( idx ) );
    if ( value != 0 )
      existingLabels.insert( value );

    if ( value > maxLabel )
      {
      maxLabel = value;
      }
    }

  // Sort out how many labels there are and allocate smallest type that can accomodate them
  TypedArray::SmartPtr resultData;
  if ( maxLabel < 256 )
    resultData = TypedArray::Create( TYPE_BYTE, nPixels );
  else if ( maxLabel < 65536 )
    resultData = TypedArray::Create( TYPE_USHORT, nPixels );
  else
    resultData = TypedArray::Create( TYPE_UINT, nPixels );
  resultData->SetDataClass( DATACLASS_LABEL );
  resultData->ClearArray();

  // Run over all labels and compound the per-label eroded label maps
  for ( std::set<unsigned int>::const_iterator it = existingLabels.begin(); it != existingLabels.end(); ++it )
    {
    // Use EDT to erode
    TypedArray::SmartPtr erodedData = UniformDistanceMap<Types::Coordinate>( volume, DistanceMap::INSIDE | DistanceMap::VALUE_EXACT, *it ).Get()->GetData();
    erodedData->Binarize( erodeBy + 0.5 );

    // Compound into final map
    for ( size_t idx = 0; idx < nPixels; ++idx )
      {
      if ( erodedData->ValueAt( idx ) > 0 )
	resultData->Set( *it, idx );
      }
    }
  
  return resultData;
}

TypedArray::SmartPtr
UniformVolumeMorphologicalOperators::GetDilatedByDistance( const Types::Coordinate dilateBy ) const
{
  if ( !this->m_UniformVolume->GetData() )
    return TypedArray::SmartPtr( NULL );
  
  TypedArray::SmartPtr dilatedData = UniformDistanceMap<Types::Coordinate>( *(this->m_UniformVolume) ).Get()->GetData();
  dilatedData->Binarize( dilateBy + 0.5 );
  dilatedData->Rescale( -1 /*scale*/, +1 /*offset*/ ); // this is binary inversion, 0->1, 1->0
  dilatedData->SetDataClass( DATACLASS_LABEL );

  return dilatedData->Convert( TYPE_BYTE );
}


} // namespace cmtk
