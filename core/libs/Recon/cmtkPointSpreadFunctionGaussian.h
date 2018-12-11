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

#ifndef __cmtkPointSpreadFunctionGaussian_h_included_
#define __cmtkPointSpreadFunctionGaussian_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup Recon */
//@{
namespace
PointSpreadFunctions
{

/// Gaussian point-spread function.
class Gaussian
{
public:
  /// This class.
  typedef Gaussian Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Constructor.
  Gaussian( const Vector3D& pixelSize )
  {
    this->m_Radius = 0.5 * pixelSize;
  }
  
  /// Get truncation radius.
  Types::Coordinate GetTruncationRadius( const int dim ) const
  {
    return 4 * this->m_Radius[dim];
  }

  /// Get the weight for a neighbor based on its radius from the kernel center.
  Types::Coordinate GetWeight( const int dim, const Types::Coordinate r ) const
  {
    Types::Coordinate rAbsRel = fabs( r / this->m_Radius[dim] );
    if ( rAbsRel <= 4 )
      {
      rAbsRel *= (2.3548 / 2.0); // go from HWHM to FWHM to sigma
      return exp( -rAbsRel*rAbsRel );
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
