/*
//
//  Copyright 2012 SRI International
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

#ifndef __cmtkFitRigidToLandmarks_h_included_
#define __cmtkFitRigidToLandmarks_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkLandmarkPairList.h>
#include <Base/cmtkMatrix3x3.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Fit affine transformation to a set of landmark pairs.
 */
class FitRigidToLandmarks
{
public:
  /// This class.
  typedef FitRigidToLandmarks Self;

  /// Constructor.
  FitRigidToLandmarks( const LandmarkPairList& landmarkPairs );

  /// Return the rigid transformation.
  AffineXform::SmartPtr GetRigidXform()
  {
    return this->m_RigidXform;
  }
  
  /// Return the constant rigid transformation.
  AffineXform::SmartConstPtr GetRigidXform() const
  {
    return this->m_RigidXform;
  }
  
private:
  /// The fitted transformation.
  AffineXform::SmartPtr m_RigidXform;
};

} // namespace

#endif // #ifndef __cmtkFitRigidToLandmarks_h_included_
