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

#ifndef __cmtkClassStreamAffineXform_h_included_
#define __cmtkClassStreamAffineXform_h_included_

#include <cmtkconfig.h>

#include <IO/cmtkClassStreamInput.h>
#include <IO/cmtkClassStreamOutput.h>

#include <Base/cmtkAffineXform.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Write affine transformation object.
ClassStreamOutput& operator << ( ClassStreamOutput& stream, const AffineXform& affineXform );

/// Read affine transformation.
ClassStreamInput& operator >> ( ClassStreamInput& stream, AffineXform::SmartPtr& affineXform );

/// Read affine transformation.
ClassStreamInput& operator >> ( ClassStreamInput& stream, AffineXform& affineXform );

//@}

} // namespace cmtk

#endif // #ifndef __cmtkClassStreamAffineXform_h_included_
