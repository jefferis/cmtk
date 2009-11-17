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
//  $Revision: 808 $
//
//  $LastChangedDate: 2009-11-17 15:07:46 -0800 (Tue, 17 Nov 2009) $
//
//  $LastChangedBy: torstenrohlfing $
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

} // namespace MathUtil

} // namespace cmtk

