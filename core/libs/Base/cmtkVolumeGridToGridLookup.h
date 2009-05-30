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

#ifndef __cmtkVolumeGridToGridLookup_h_included_
#define __cmtkVolumeGridToGridLookup_h_included_

#include <cmtkconfig.h>

#include <cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Class for grid to grid lookup tables.
 * This class is only for internal use by the Resmple function(s) of
 * the UniformVolume class. The lookup table basically records for
 * each pixel in the target grid, which pixels in the source grid it 
 * depends on, and with what relative weight each source pixel contributes
 * to the target pixel. The contributions are computed as relative overlaps
 * of boxcar-shaped pixel profiles.
 */
class VolumeGridToGridLookup
{
public:
  /// Constructor: takes original and new image grids.
  VolumeGridToGridLookup( const UniformVolume& fromGrid, const UniformVolume& toGrid );

  /// Get number of source pixels that contribute to the given target pixel.
  int GetSourceCount( const int dim, const int idx ) const
  {
    return this->m_SourceCount[dim][idx];
  }

  /// Get index of first source pixel that contributes to the given target pixel.
  int GetFromIndex( const int dim, const int idx ) const
  {
    return this->m_FromIndex[dim][idx];
  }

  /// Get weight with which a given source pixel contributes to the given target pixel.
  Types::Coordinate GetWeight( const int dim, const int idx, const int fromIdx ) const
  {
    return this->m_Weight[dim][idx][fromIdx];
  }

  /// Length (width) of a given target pixel.
  Types::Coordinate GetLength( const int dim, const int idx ) const
  {
    return this->m_Length[dim][idx];
  }

private:
  /// Array of arrays of numbers of source pixels that contribute to the given target pixel.
  std::vector< std::vector< int > > m_SourceCount;

  /// Array of arrays of first source pixels that contributes to the given target pixels.
  std::vector< std::vector< int > > m_FromIndex;

  /// Array of arrays of weight arrays.
  std::vector< std::vector< std::vector<Types::Coordinate> > > m_Weight;

  /// Array of array of target pixel lengths.
  std::vector< std::vector< Types::Coordinate > > m_Length;
};

//@}

} // namespace cmtk

#endif // #define __cmtkUniformVolume_h_included_
