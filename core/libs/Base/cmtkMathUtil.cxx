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

#include <cmtkMathUtil.h>

#include <cmtkTypes.h>
#include <cmtkMatrix.h>

#include <assert.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

namespace
MathUtil
{

const FPInitializeUnion FPInitializeNaN = {
#if WORDS_BIGENDIAN
  { 0x7fffffff, 0xffffffff }
#else
  { 0xffffffff, 0x7fffffff }
#endif
};

const FPInitializeUnion FPInitializeInf = {
#if WORDS_BIGENDIAN
  { 0x7f800000, 0x00000000 }
#else
  { 0x00000000, 0x7f800000 }
#endif
};

const void* FPInitializeNaN_P = &FPInitializeNaN;
const void* FPInitializeInf_P = &FPInitializeInf;

#if WORDS_BIGENDIAN
const void* FPInitializeNaN_fP = &FPInitializeNaN.f[0];
const void* FPInitializeInf_fP = &FPInitializeInf.f[0];
#else
const void* FPInitializeNaN_fP = &FPInitializeNaN.f[1];
const void* FPInitializeInf_fP = &FPInitializeInf.f[1];
#endif


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

