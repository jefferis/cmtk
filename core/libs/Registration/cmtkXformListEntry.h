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

#ifndef __cmtkXformListEntry_h_included_
#define __cmtkXformListEntry_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkXform.h"
#include "Base/cmtkAffineXform.h"
#include "Base/cmtkWarpXform.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{
/// An entry in a transformation sequence.
class XformListEntry 
{
public:
  /// Constructor.
  XformListEntry( const Xform::SmartPtr& xform = Xform::SmartPtr::Null, const bool inverse = false, const Types::Coordinate globalScale = 1.0 );
  
  /// Destructor.
  ~XformListEntry();
  
  /// The actual transformation.
  const Xform::SmartPtr m_Xform;
  
  /// The actual inverse if transformation is affine.
  const AffineXform* InverseAffineXform;
  
  /// The actual transformation as spline warp.
  const WarpXform* m_WarpXform;
  
  /// Apply forward (false) or inverse (true) transformation.
  bool Inverse;
  
  /// Global scale for normalizing the Jacobian.
  Types::Coordinate GlobalScale;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkXformListEntry_h_included_
