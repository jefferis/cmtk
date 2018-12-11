/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkMathFunctionWrappers_h_included_
#define __cmtkMathFunctionWrappers_h_included_

#include <cmtkconfig.h>

#include <math.h>

namespace
cmtk
{

/** Wrappers for math.h functions.
 * These are necessary because on some platforms (e.g., Solaris), the
 * math functions are declared extern "C"; so until we find a better way
 * to address this, we wrap these functions.
 */

namespace
Wrappers
{

/// Log function.
double Log( const double x );

/// Exponential function.
double Exp( const double x );

/// Square root function.
double Sqrt( const double x );

/// Abs function.
double Abs( const double x );

/// Trunc function.
double Trunc( const double x );

/// Square function.
double Square( const double x );

/// Logit function.
double Logit( const double x );

/// Logistic function.
double Logistic( const double x );

} // namespace Wrappers

} // namespace cmtk

#endif // #ifndef __cmtkMathFunctionWrappers_h_included_
