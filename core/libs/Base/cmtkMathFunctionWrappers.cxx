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

#include "cmtkMathFunctionWrappers.h"

namespace
cmtk
{

namespace
Wrappers
{

double
Log( const double x )
{
  return log( x );
}

double
Exp( const double x )
{
  return exp( x );
}

double
Sqrt( const double x )
{
  return sqrt( x );
}

double
Abs( const double x )
{
  return fabs( x );
}

double
Trunc( const double x )
{
#ifdef _MSC_VER
  return static_cast<double>( static_cast<long int>( x ) );
#else
  return trunc( x );
#endif
}

double Square( const double x )
{
  return x*x; 
}

double Logit( const double x )
{
  return log(x / (1.0-x)); 
}

double Logistic( const double x )
{
  return 1.0/(1.0+exp(-x));
}

} // namespace Wrappers

} // namespace cmtk

