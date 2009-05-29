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

#ifndef __cmtkDirectionSet_h_included_
#define __cmtkDirectionSet_h_included_

#include <cmtkconfig.h>

#include <cmtkAffineXform.h>
#include <cmtkVector.h>
#include <cmtkSmartPtr.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// A set of directions in n-dimensional space.
class DirectionSet :
  /// This is a vector of coordinate vectors.
  public std::vector<CoordinateVector::SmartPtr>
{
public:
  /// Smart pointer to DirectionSet.
  typedef SmartPointer<DirectionSet> SmartPtr;

  /// Get dimension of direction space.
  unsigned int GetDimension() const { return Dimension; }

  /// Get number of directions.
  unsigned int GetNumberOfDirections() const { return this->size(); }

  /** Normalizes each direction vector to have max norm = value.
   */
  void NormalizeMaxNorm( const double value = 1.0 );

  /** Normalizes each direction vector to have euclid norm = value.
   */
  void NormalizeEuclidNorm( const double value = 1.0 );

  /** Reorients the direction vectors based on an affine transformation.
   */
  void Reorient( const AffineXform* affineXform );

  /// Default constructor.
  DirectionSet() { Dimension = 0; }

  /// Allocation constructor.
  DirectionSet( const unsigned int dimension ) { Dimension = dimension; }

  /// Destructor.
  ~DirectionSet() {}

private:
  /// Dimension of direction space.
  unsigned int Dimension;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDirectionSet_h_included_
