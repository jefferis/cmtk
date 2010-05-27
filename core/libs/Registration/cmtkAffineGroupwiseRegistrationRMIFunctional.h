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

#ifndef __cmtkAffineGroupwiseRegistrationRMIFunctional_h_included_
#define __cmtkAffineGroupwiseRegistrationRMIFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkGroupwiseRegistrationRMIFunctional.h>

#include <cmtkSmartPtr.h>
#include <cmtkAffineXform.h>
#include <cmtkClassStream.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Functional for affine groupwise registration using "RMI" metric..
 */
class AffineGroupwiseRegistrationRMIFunctional : 
  /// Inherit general RMI-based groupwise registration.
  public GroupwiseRegistrationRMIFunctional<AffineXform>
{
public:
  /// Type of parent class.
  typedef GroupwiseRegistrationRMIFunctional<AffineXform> Superclass;

  /// Type of this class.
  typedef AffineGroupwiseRegistrationRMIFunctional Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  friend ClassStream& operator<<( ClassStream& stream, const AffineGroupwiseRegistrationRMIFunctional& func );
  friend ClassStream& operator>>( ClassStream& stream, AffineGroupwiseRegistrationRMIFunctional& func );
};

/// Class stream write function.
ClassStream& operator<<( ClassStream& stream, const AffineGroupwiseRegistrationRMIFunctional& func );

/// Class stream read function.
ClassStream& operator>>( ClassStream& stream, AffineGroupwiseRegistrationRMIFunctional& func );

//@}

} // namespace cmtk

#endif // #ifndef __cmtkAffineGroupwiseRegistrationRMIFunctional_h_included_
