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

#include "cmtkMathUtil.h"

#include <Base/cmtkTypes.h>
#include <System/cmtkConsole.h>
#include "Numerics/ibetaf.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{

double
MathUtil::TStatFromCorrelation
( const double r, const size_t df )
{
  return r * sqrt( df / (1-r*r) ); 
}

double 
MathUtil::ProbabilityFromTStat
( const double t, const size_t df )
{
  double stat;
  if ( df == 0.0 )
    stat = 0.0;
  else if ( t == 0.0 )
    stat = 1.0;
  else
    stat = df/(df+t*t);

  return alglib::incompletebeta( 0.5*df, 0.5, stat );
}

} // namespace cmtk
