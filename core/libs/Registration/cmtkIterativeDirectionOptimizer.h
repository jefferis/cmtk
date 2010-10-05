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

#ifndef __cmtkIterativeDirectionOptimizer_h_included_
#define __cmtkIterativeDirectionOptimizer_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkOptimizer.h>

namespace 
cmtk
{

/** \addtogroup Registration */
//@{

/** Iterative direction optimizer.
 * This class implements a simple search algorithm that will iterate through
 * all dimensions of the parameter space and search for the maximum of the
 * target function along each dimension. If no further improvements are
 * possible, the step size is reduced and the search continues until a given
 * minimum step size is reached.
 */
class IterativeDirectionOptimizer : 
  /// Inherit generic optimizer features.
  public Optimizer 
{
public:
  /** Constructor.
   * Hand functional and callback to parent class and initialize local 
   * variables.
   *@param stepFactor Factor by which the search step size is decreased.
   */
  IterativeDirectionOptimizer ( const Types::Coordinate stepFactor = 0.5 )
  { StepFactor = stepFactor; };
 
  /** Perform the optimization.
   */
  virtual CallbackResult Optimize( CoordinateVector&, const Types::Coordinate = 1, const Types::Coordinate = 0 );

private:
  /** Search step factor.
   * This variable determines the factor by which to decrease the search step
   * size if no further improvement is possible at a certain resolution.
   * Reasonable values are in the range 0 < StepFactor < 1. For most cases,
   * a value of 0.5 has been found to provide optimum accuracy and performance.
   */
  Types::Coordinate StepFactor;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkIterativeDirectionOptimizer_h_included_
