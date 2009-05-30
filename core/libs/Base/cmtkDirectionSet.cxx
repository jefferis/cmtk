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

#include <cmtkDirectionSet.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void 
DirectionSet::NormalizeMaxNorm( const double value )
{
  // For Each Direction
  for( size_t index = 0; index < this->GetNumberOfDirections(); index++ )
    {
    CoordinateVector::SmartPtr direction = (*this)[index];
    (*direction) *= value / direction->MaxNorm();
    }
}

void 
DirectionSet::NormalizeEuclidNorm(const double value)
{
  // For Each Direction
  for( size_t index = 0; index < this->GetNumberOfDirections(); index++)
    {
    CoordinateVector::SmartPtr direction = (*this)[index];
    
    // Find Euclidean Norm
    Types::Coordinate euclid = direction->EuclidNorm();
    
    // Divide Every Component by the Euclidean Norm
    // and multiply by normalization value..
    *direction *= value / euclid;
    }
}

void
DirectionSet::Reorient( const AffineXform* affineXform )
{
  // For Each Direction
  for( size_t index = 0; index < this->GetNumberOfDirections(); index++ )
    {
    CoordinateVector::SmartPtr direction = (*this)[index];
    
    // For each point
    Types::Coordinate *point = direction->Elements;
    for( size_t i = 0; i < this->GetDimension(); i += 3, point += 3)
      {
      // Transform Point
      affineXform->RotateScaleShear( point, point );
      }
    }
}

} // namespace cmtk
