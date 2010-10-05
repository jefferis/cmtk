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

#ifndef __cmtkClassStreamMultiChannelRegistration_h_included_
#define __cmtkClassStreamMultiChannelRegistration_h_included_

#include <IO/cmtkClassStream.h>

#include <Registration/cmtkAffineMultiChannelRegistrationFunctional.h>
#include <Registration/cmtkSplineWarpMultiChannelRegistrationFunctional.h>

namespace 
cmtk
{

/** \addtogroup IO */
//@{

/** Write file names and transformations from affine multi-channel registration functional. */
template<class TMetricFunctionalType>
ClassStream& operator << ( ClassStream& stream, const AffineMultiChannelRegistrationFunctional<TMetricFunctionalType>& functional );

/** Read file names and transformations from archive to multi-channel affine registration functional. */
template<class TMetricFunctionalType>
ClassStream& operator >> ( ClassStream& stream, AffineMultiChannelRegistrationFunctional<TMetricFunctionalType>& functional );

/** Write file names and transformations from spline warp multi-channel registration functional. */
template<class TMetricFunctionalType>
ClassStream& operator << ( ClassStream& stream, const SplineWarpMultiChannelRegistrationFunctional<TMetricFunctionalType>& functional );

//@}

} // namespace cmtk

#include "cmtkClassStreamMultiChannelRegistration.txx"

#endif // #ifndef __cmtkClassStreamMultiChannelRegistration_h_included_
