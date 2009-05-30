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

#ifndef __cmtkPointSpreadFunctionBox_h_included_
#define __cmtkPointSpreadFunctionBox_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup Unstable */
//@{

/// Point spread functions for iterative deblurring.
namespace
PointSpreadFunctions
{

/// Box point-spread function.
class Box
{
public:
  /// This class.
  typedef Box Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Constructor.
  Box( const Vector3D& pixelSize )
  {
    this->m_Radius = 0.5 * pixelSize;
  }
  
  /// Get truncation radius.
  Types::Coordinate GetTruncationRadius( const int dim ) const
  {
    return this->m_Radius.XYZ[dim];
  }

  /// Get the weight for a neighbor based on its radius from the kernel center.
  Types::Coordinate GetWeight( const int dim, const Types::Coordinate r ) const
  {
    if ( fabs( r ) <= this->m_Radius.XYZ[dim] )
      {
      return 1.0;
      }
    return 0.0;
  }
  
private:
  /// Kernel radius.
  Vector3D m_Radius;
};

} // namespace PointSpreadFunctions

//@}

} // namespace cmtk

#endif
