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

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class TAccumulator>
ScalarImage*
DataGrid::ComputeProjection
( const int axis ) const
{
  unsigned int dims[2], depth, offset, incX, incY, incZ;

  switch ( axis ) 
    {
    case AXIS_X:
      dims[0] = this->m_Dims[1];
      dims[1] = this->m_Dims[2];
      depth = this->m_Dims[0];
      offset = 0;
      incX = this->m_Dims[0];
      incY = this->m_Dims[0] * this->m_Dims[1];
      incZ = 1;
      break;
    case AXIS_Y:
      dims[0] = this->m_Dims[0];
      dims[1] = this->m_Dims[2];
      depth = this->m_Dims[1];
      offset = 0;
      incX = 1;
      incY = this->m_Dims[0] * this->m_Dims[1];
      incZ = this->m_Dims[0];
      break;
    case AXIS_Z:
    default:
      dims[0] = this->m_Dims[0];
      dims[1] = this->m_Dims[1];
      depth = this->m_Dims[2];
      offset = 0;
      incX = 1;
      incY = this->m_Dims[0];
      incZ = this->m_Dims[0] * this->m_Dims[1];
      break;
    }
  
  const TypedArray* data = this->GetData();
  TypedArray::SmartPtr projectData = TypedArray::Create( data->GetType(), dims[0] * dims[1] );
  
  for ( unsigned int y = 0; y < dims[1]; ++y ) 
    {
    size_t projectOffset = y * dims[0];
    unsigned int offsetY = offset + incY;
    for ( unsigned int x = 0; x < dims[0]; ++x, ++projectOffset ) 
      {
      unsigned int offsetX = offset + incX;
      
      TAccumulator accu;
      Types::DataItem item;
      for ( unsigned int z = 0; z < depth; ++z, offset += incZ ) 
	{
       	if ( data->Get( item, offset ) ) 
	  {
	  accu.AddValue( item );
	  }
	}
      projectData->Set( accu.GetResult(), projectOffset );
      
      offset = offsetX;
      }
    offset = offsetY;
    }
  
  ScalarImage* projectImage = new ScalarImage( dims[0], dims[1] );
  projectImage->SetPixelData( TypedArray::SmartPtr( projectData ) );
  
  return projectImage;
}

template<class TData>
TData
DataGrid::TrilinearInterpolation
( const TData* dataPtr, 
  const int X, const int Y, const int Z,
  const Self::SpaceVectorType& Location, const Types::Coordinate* from, 
  const Types::Coordinate* to ) const
{
  const size_t offset = X+this->m_Dims[0]*(Y+this->m_Dims[1]*Z);
  const TData* data = dataPtr + offset;

  const Types::Coordinate deltaX=1.0/(to[0]-from[0]), deltaY=1.0/(to[1]-from[1]), deltaZ=1.0/(to[2]-from[2]);
  
  const Types::Coordinate revX = deltaX*(Location[0]-from[0]);
  const Types::Coordinate revY = deltaY*(Location[1]-from[1]);
  const Types::Coordinate revZ = deltaZ*(Location[2]-from[2]);
  const Types::Coordinate offsX = 1-revX;
  const Types::Coordinate offsY = 1-revY;
  const Types::Coordinate offsZ = 1-revZ;
  
  return static_cast<TData>( offsZ*(offsY*(offsX*data[0] + revX*data[nextI])+ 
				    revY*(offsX*data[nextJ]+ revX*data[nextIJ]))+
			     revZ*(offsY*(offsX*data[nextK]+ revX*data[nextIK])+ 
				   revY*(offsX*data[nextJK]+ revX*data[nextIJK])));
}

template<class TData,class TOutputIterator>
inline void
DataGrid
::TrilinearInterpolation
( TOutputIterator result, const std::vector<TData*>& dataPtr, const int x, const int y, const int z,
  const Types::Coordinate fracX, const Types::Coordinate fracY, const Types::Coordinate fracZ ) const
{
  const size_t offset = x + this->m_Dims[0] * ( y + this->m_Dims[1] * z);

  const Types::Coordinate offsX = 1.0-fracX;
  const Types::Coordinate offsY = 1.0-fracY;
  const Types::Coordinate offsZ = 1.0-fracZ;

  const Types::Coordinate tmp0 = offsZ * offsY;
  const Types::Coordinate w0 = tmp0 * offsX;
  const Types::Coordinate w1 = tmp0 * fracX;

  const Types::Coordinate tmp1 = offsZ * fracY;
  const Types::Coordinate w2 = tmp1 * offsX;
  const Types::Coordinate w3 = tmp1 *  fracX;

  const Types::Coordinate tmp2 = fracZ * offsY;
  const Types::Coordinate w4 = tmp2 * offsX;
  const Types::Coordinate w5 = tmp2 * fracX;

  const Types::Coordinate tmp3 = fracZ *  fracY;
  const Types::Coordinate w6 = tmp3 * offsX;
  const Types::Coordinate w7 = tmp3 * fracX;
  
  for ( size_t i = 0; i < dataPtr.size(); ++i, ++result )
    {    
    const TData* data = dataPtr[i] + offset;
    *result = static_cast<TData>( w0 * data[0] +     w1 * data[nextI] +  w2 * data[nextJ] +  w3 * data[nextIJ] +
				  w4 * data[nextK]+  w5 * data[nextIK] + w6 * data[nextJK] + w7 * data[nextIJK] );
    }
}

} // namespace cmtk
