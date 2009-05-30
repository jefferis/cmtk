/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkReformatVolume.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

TypedArray* 
ReformatVolume::ReformatPushForward
( const UniformVolume* floating, XformList& targetToRef, const UniformVolume* reference )
{
  const int* dims = floating->GetDims();
  int xyzCP[3];

  TypedArray* result = TypedArray::Create( floating->GetData()->GetType(), reference->GetNumberOfPixels() );
  result->Fill( 0.0 );

  for ( int z = 0; z < dims[2]; ++z )
    for ( int y = 0; y < dims[1]; ++y )
      for ( int x = 0; x < dims[0]; ++x )
	{
	Vector3D v( floating->GetGridLocation( x, y, z ) );
	targetToRef.Apply( v );
	if ( reference->GetClosestGridPointIndex( v, xyzCP ) )
	  {
	  Types::DataItem value;
	  if ( floating->GetDataAt( value, x, y, z ) )
	    {
	    const size_t offset = reference->GetOffsetFromIndex( xyzCP[0], xyzCP[1], xyzCP[2] );
	    result->Set( value, offset );
	    }
	  }
	
	}
  return result;
}

TypedArray* 
ReformatVolume::ReformatPushForwardAccumulate
( const UniformVolume* floating, XformList& targetToRef, const UniformVolume* reference )
{
  const int* dims = floating->GetDims();
  int xyzCP[3];

  TypedArray* result = TypedArray::Create( floating->GetData()->GetType(), reference->GetNumberOfPixels() );
  result->Fill( 0.0 );

  for ( int z = 0; z < dims[2]; ++z )
    for ( int y = 0; y < dims[1]; ++y )
      for ( int x = 0; x < dims[0]; ++x )
	{
	Vector3D v( floating->GetGridLocation( x, y, z ) );
	targetToRef.Apply( v );
	if ( reference->GetClosestGridPointIndex( v, xyzCP ) )
	  {
	  Types::DataItem value;
	  if ( floating->GetDataAt( value, x, y, z ) )
	    {
	    const size_t offset = reference->GetOffsetFromIndex( xyzCP[0], xyzCP[1], xyzCP[2] );
            Types::DataItem valuePrev;
	    result->Get( valuePrev, offset );
	    result->Set( std::max( value, valuePrev ), offset );
	    }
	  }
	}
  
  return result;
}

} // namespace cmtk
