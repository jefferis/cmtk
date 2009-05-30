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

#include <cmtkUniformVolume.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

TypedArray*
UniformVolume::GetDataGaussFiltered( const Types::Coordinate stdDev ) const
{
  const Types::Coordinate stdDevPixelX = stdDev / this->Delta[0];
  const Types::Coordinate stdDevPixelY = stdDev / this->Delta[1];
  const Types::Coordinate stdDevPixelZ = stdDev / this->Delta[2];

  const unsigned int stdDevDiscreteX = static_cast<unsigned int>( ceil( stdDevPixelX ) );
  const unsigned int stdDevDiscreteY = static_cast<unsigned int>( ceil( stdDevPixelY ) );
  const unsigned int stdDevDiscreteZ = static_cast<unsigned int>( ceil( stdDevPixelZ ) );

  const unsigned int filterLengthX = std::min<unsigned int>( this->Dims[0], 3 * stdDevDiscreteX + 1 );
  const unsigned int filterLengthY = std::min<unsigned int>( this->Dims[1], 3 * stdDevDiscreteY + 1 );
  const unsigned int filterLengthZ = std::min<unsigned int>( this->Dims[2], 3 * stdDevDiscreteZ + 1 );

  std::vector<Types::DataItem> filterX( filterLengthX );
  for ( unsigned int x=0; x < filterLengthX; ++x ) 
    {
    filterX[x] = 1.0/(sqrt(2*M_PI) * stdDevPixelX) * exp( -MathUtil::Square( 1.0 * x / stdDevPixelX ) / 2 );
    }
  
  std::vector<Types::DataItem> filterY( filterLengthY );
  for ( unsigned int y=0; y < filterLengthY; ++y ) 
    {
    filterY[y] = 1.0/(sqrt(2*M_PI) * stdDevPixelY) * exp( -MathUtil::Square( 1.0 * y / stdDevPixelY ) / 2);
    }
  
  std::vector<Types::DataItem> filterZ( filterLengthZ );
  for ( unsigned int z=0; z < filterLengthZ; ++z ) 
    {
    filterZ[z] = 1.0/(sqrt(2*M_PI) * stdDevPixelZ) * exp( -MathUtil::Square( 1.0 * z / stdDevPixelZ ) / 2);
    }
  
  TypedArray *result = this->GetFilteredData( filterX, filterY, filterZ );
  
  return result;
}

} // namespace cmtk
