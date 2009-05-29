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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkVolume.h>
#include <cmtkVolumeGridToGridLookup.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

VolumeGridToGridLookup
::VolumeGridToGridLookup( const UniformVolume& fromVolume, const UniformVolume& toVolume )
  : m_SourceCount( 3 ), m_FromIndex( 3 ), m_Weight( 3 ), m_Length( 3 )
{
  for ( int dim = 0; dim < 3; ++dim )
    {
    const Types::Coordinate fromGridDelta = fromVolume.Delta[dim];
    const Types::Coordinate toGridDelta = toVolume.Delta[dim];

    this->m_SourceCount[dim].resize( toVolume.Dims[dim]+1 );
    this->m_FromIndex[dim].resize( toVolume.Dims[dim]+1 );
    this->m_Weight[dim].resize( toVolume.Dims[dim]+1 );
    this->m_Length[dim].resize( toVolume.Dims[dim]+1 );

    std::vector<Types::Coordinate> weightList( fromVolume.Dims[dim] );
    
    int fromIdx = 0;
    for ( int toIdx = 0; toIdx < toVolume.Dims[dim]; ++toIdx ) 
      {
      const Types::Coordinate toGridLo = std::max<Types::Coordinate>( 0.0, (toIdx-0.5) * toGridDelta );
      const Types::Coordinate toGridHi = std::min<Types::Coordinate>( toVolume.Size[dim], (0.5+toIdx) * toGridDelta );
      this->m_Length[dim][toIdx] = toGridHi - toGridLo;

      Types::Coordinate fromGridHi = std::min<Types::Coordinate>( toVolume.Size[dim], (0.5+fromIdx) * fromGridDelta );
      while ( toGridLo>=fromGridHi )
	{
	++fromIdx;
	fromGridHi += fromGridDelta;
	}
      this->m_FromIndex[dim][toIdx] = fromIdx;
      fromGridHi = std::min<Types::Coordinate>( fromVolume.Size[dim], fromGridHi );
      
      int idx = 0;
      Types::Coordinate fromGridLo = std::max<Types::Coordinate>( 0.0, (fromIdx-0.5) * fromGridDelta );
      for ( int p = fromIdx; (p < fromVolume.Dims[dim]) && (fromGridLo < toGridHi); ++p, ++idx )
	{
	weightList[idx] = MathUtil::Intersect( toGridLo, toGridHi, fromGridLo, fromGridHi );
	fromGridLo = (p+0.5) * fromGridDelta;	
	fromGridHi += fromGridDelta;
	}
      
      this->m_SourceCount[dim][toIdx] = idx;
      this->m_Weight[dim][toIdx].resize( idx );
      for ( int i = 0; i < idx; ++i )
	this->m_Weight[dim][toIdx][i] = weightList[i];
      }
    
    this->m_Weight[dim][toVolume.Dims[dim]].resize(0);
    }
}

} // namespace cmtk
