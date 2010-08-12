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

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TInterpolator, class Fct> 
TypedArray::SmartPtr
ReformatVolume::Reformat
( const UniformVolume* target, cmtk::XformList& targetToRef, const UniformVolume* reference, cmtk::XformList& refToFloat, Fct& fct, TInterpolator& interpolator )
{
  const DataGrid::IndexType& dims = target->GetDims();

  TypedArray::SmartPtr result = TypedArray::Create( fct.GetDataType( reference, interpolator ), target->GetNumberOfPixels() );
  if ( fct.UsePaddingValue )
    result->SetPaddingValue( fct.PaddingValue );
  const TypedArray* targetData = target->GetData();
  
  Progress::Begin( 0, dims[2], 1, "Volume reformatting" );
  
#pragma omp parallel for
  for ( int z = 0; z < dims[2]; z++ ) 
    {
    Vector3D vRef;
    Types::DataItem value, targetValue;
    size_t offset = z * dims[0] * dims[1];
    Progress::SetProgress( z );
    
    for ( int y = 0; y < dims[1]; y++ ) 
      {
      for ( int x = 0; x < dims[0]; x++, offset++ ) 
	{
	if ( !targetData || (targetData->Get( targetValue, offset ) && (targetValue != 0))) 
	  {
	  vRef = target->GetGridLocation( x, y, z );
	  if ( targetToRef.ApplyInPlace( vRef ) && fct( value, vRef, refToFloat, interpolator ) ) 
	    {
	    result->Set( value, offset );
	    } 
	  else
	    {
	    result->SetPaddingAt( offset );
	    }
	  } 
	else
	  {
	  result->SetPaddingAt( offset );
	  }
	}
      }
    }
  
  Progress::Done();
  return result;
}

} // namespace cmtk
