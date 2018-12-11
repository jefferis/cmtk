/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#ifndef __cmtkCompatibilityMatrix4x4_h_included_
#define __cmtkCompatibilityMatrix4x4_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkMatrix4x4.h>
#include <Base/cmtkVector.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Compatibility class for homogeneous 4x4 transformation matrix.
 * The sole purpose of this class is to take a parameter vector describing the degrees of freedom of a 
 * cmtk::Matrix4x4 object constructed prior to CMTK release 2.4. Older releases of CMTK generated matrices
 * in which scale and shear cooefficients were not fully independent, making it impossible to recover the
 * exact parameters from the matrix.
 *\see https://www.nitrc.org/tracker/index.php?func=detail&aid=7179&group_id=212&atid=877
 */
template<class T=Types::Coordinate>
class CompatibilityMatrix4x4 :
    public Matrix4x4<T>
{
public:
  /// This class.
  typedef CompatibilityMatrix4x4<T> Self;

  /// Parent class..
  typedef Matrix4x4<T> Superclass;

  /// Constructor: create matrix from parameter vector as it would have been done prior to CMTK 2.4.
  CompatibilityMatrix4x4( const CoordinateVector& dofs, const bool logScaleFactors = false );
};

} // namespace cmtk


#endif // #ifndef __cmtkCompatibilityMatrix4x4_h_included_
