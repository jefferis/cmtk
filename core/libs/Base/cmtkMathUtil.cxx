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

#include <cmtkMathUtil.h>

#include <math.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

namespace
MathUtil
{

#ifdef HAVE_FABSF
inline float absF ( const float a ) { return fabsf(a); }
#else
inline float absF ( const float a ) { return static_cast<float>( fabs(a) ); }
#endif
inline double absF ( const double a ) { return fabs(a); }

#ifdef HAVE_FMODF
inline float modF ( const float a, const float b ) { return fmodf( a, b ); }
#else
inline float modF ( const float a, const float b ) 
{ return static_cast<float>( fmod( a, b ) ); }
#endif

inline double modF ( const double a, const double b ) { return fmod( a, b ); }

} // namespace MathUtil

} // namespace cmtk

