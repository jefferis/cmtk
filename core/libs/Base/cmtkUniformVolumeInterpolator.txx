/*
//
//  Copyright 2016 Google, Inc.
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

namespace
cmtk
{

/** \addtogroup Base */
//@{
template <class TInterpolationFunction>
bool
UniformVolumeInterpolator<TInterpolationFunction>
::GetDataAt(const Vector3D& v, Types::DataItem& value) const
{
  Types::Coordinate lScaled[3];
  Types::GridIndexType imageGridPoint[3];
  for ( int n = 0; n < 3; ++n )
    {
    lScaled[n] = (v[n] - this->m_VolumeOffset[n]) / this->m_VolumeDeltas[n];
    imageGridPoint[n] = (int) floor( lScaled[n] );
    if ( ( imageGridPoint[n] < 0 ) || ( imageGridPoint[n] >= this->m_VolumeDims[n]-1 ) )
      return false;
    }

  const Types::GridIndexType xx = imageGridPoint[0] + 1 - TInterpolationFunction::RegionSizeLeftRight;
  const Types::GridIndexType yy = imageGridPoint[1] + 1 - TInterpolationFunction::RegionSizeLeftRight;
  const Types::GridIndexType zz = imageGridPoint[2] + 1 - TInterpolationFunction::RegionSizeLeftRight;

  Types::Coordinate interpolationWeights[3][2 * TInterpolationFunction::RegionSizeLeftRight];

  for ( int n = 0; n < 3; ++n )
    {
    const Types::Coordinate relative = lScaled[n] - imageGridPoint[n];
    for ( Types::GridIndexType m = 1-TInterpolationFunction::RegionSizeLeftRight; m <= TInterpolationFunction::RegionSizeLeftRight; ++m )
      {
      interpolationWeights[n][m+TInterpolationFunction::RegionSizeLeftRight-1] = TInterpolationFunction::GetWeight( m, relative );
      }
    }

  const Types::GridIndexType iMin = std::max<Types::GridIndexType>( 0, -xx );
  const Types::GridIndexType iMax = std::min<Types::GridIndexType>( 2 * TInterpolationFunction::RegionSizeLeftRight, this->m_VolumeDims[0] - xx );

  const Types::GridIndexType jMin = std::max<Types::GridIndexType>( 0, -yy );
  const Types::GridIndexType jMax = std::min<Types::GridIndexType>( 2 * TInterpolationFunction::RegionSizeLeftRight, this->m_VolumeDims[1] - yy );

  const Types::GridIndexType kMin = std::max<Types::GridIndexType>( 0, -zz );
  const Types::GridIndexType kMax = std::min<Types::GridIndexType>( 2 * TInterpolationFunction::RegionSizeLeftRight, this->m_VolumeDims[2] - zz );

  Types::DataItem interpolatedData = 0;
  Types::Coordinate totalWeight = 0;

  for ( Types::GridIndexType k = kMin; k < kMax; ++k )
    {
    for ( Types::GridIndexType j = jMin; j < jMax; ++j )
      {
      const Types::Coordinate weightJK = interpolationWeights[1][j] * interpolationWeights[2][k];
      size_t offset = this->GetOffsetFromIndex( xx + iMin, yy + j, zz + k );
      for ( Types::GridIndexType i = iMin; i < iMax; ++i, ++offset )
        {
	const Types::DataItem data = this->m_VolumeDataArray[offset];
	if ( finite( data ) )
          {
	  const Types::Coordinate weightIJK = interpolationWeights[0][i] * weightJK;
	  interpolatedData += static_cast<Types::DataItem>( data * weightIJK );
          totalWeight += weightIJK;
          }
        }
      }
    }

  if ( totalWeight == 0 )
    return false;
  else
    value = static_cast<Types::DataItem>( interpolatedData / totalWeight );

  return true;
}

template <class TInterpolationFunction>
Types::DataItem
UniformVolumeInterpolator<TInterpolationFunction>
::GetDataDirect( const Types::GridIndexType* imageGridPoint, const Types::Coordinate* insidePixel ) const
{
  Types::Coordinate interpolationWeights[3][2 * TInterpolationFunction::RegionSizeLeftRight];

  for ( int n = 0; n < 3; ++n )
    {
    for ( Types::GridIndexType m = 1-TInterpolationFunction::RegionSizeLeftRight; m <= TInterpolationFunction::RegionSizeLeftRight; ++m )
      {
      interpolationWeights[n][m+TInterpolationFunction::RegionSizeLeftRight-1] = TInterpolationFunction::GetWeight( m, insidePixel[n] );
      }
    }

  const Types::GridIndexType xx = imageGridPoint[0] + 1 - TInterpolationFunction::RegionSizeLeftRight;
  const Types::GridIndexType yy = imageGridPoint[1] + 1 - TInterpolationFunction::RegionSizeLeftRight;
  const Types::GridIndexType zz = imageGridPoint[2] + 1 - TInterpolationFunction::RegionSizeLeftRight;

  const Types::GridIndexType iMin = std::max<Types::GridIndexType>( 0, -xx );
  const Types::GridIndexType iMax = std::min<Types::GridIndexType>( 2 * TInterpolationFunction::RegionSizeLeftRight, this->m_VolumeDims[0] - xx );

  const Types::GridIndexType jMin = std::max<Types::GridIndexType>( 0, -yy );
  const Types::GridIndexType jMax = std::min<Types::GridIndexType>( 2 * TInterpolationFunction::RegionSizeLeftRight, this->m_VolumeDims[1] - yy );

  const Types::GridIndexType kMin = std::max<Types::GridIndexType>( 0, -zz );
  const Types::GridIndexType kMax = std::min<Types::GridIndexType>( 2 * TInterpolationFunction::RegionSizeLeftRight, this->m_VolumeDims[2] - zz );

  Types::DataItem interpolatedData = 0;
  Types::Coordinate totalWeight = 0;

  for ( Types::GridIndexType k = kMin; k < kMax; ++k )
    {
    for ( Types::GridIndexType j = jMin; j < jMax; ++j )
      {
      const Types::Coordinate weightJK = interpolationWeights[1][j] * interpolationWeights[2][k];
      size_t offset = (xx + iMin) + (yy + j) * this->m_NextJ + (zz + k) * this->m_NextK;
      for ( Types::GridIndexType i = iMin; i < iMax; ++i, ++offset )
        {
        const Types::DataItem data = this->m_VolumeDataArray[offset];
	if ( finite( data ) )
          {
	  const Types::Coordinate weightIJK = interpolationWeights[0][i] * weightJK;
          interpolatedData += static_cast<Types::DataItem>( data * weightIJK );
          totalWeight += weightIJK;
          }
        }
      }
    }
  
  if ( totalWeight == 0 )
    return 0;
  else
    return static_cast<Types::DataItem>( interpolatedData / totalWeight );
}

} // namespace cmtk
