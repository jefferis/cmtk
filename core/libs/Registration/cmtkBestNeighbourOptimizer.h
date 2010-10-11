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

#ifndef __cmtkBestNeighbourOptimizer_h_included_
#define __cmtkBestNeighbourOptimizer_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkOptimizer.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Best-neighbour-search optimizer.
 * This class implements a search technique introduced by Studholme et al.
 * By modifying each parameter of the search space by a certain step up- and
 * downwards, all "neighbours" of the current parameter vector are visited. For
 * each of these, the target function (functional) is evaluated. The search 
 * then continues from the parameter vector producing the maximum value. If no
 * further improvement is possible, the step size is decreased by a given 
 * factor until it reaches a lower bound.
 */
class BestNeighbourOptimizer : 
  /// Inherit generic optimizer features.
  public Optimizer 
{
public:
  /// This class.
  typedef BestNeighbourOptimizer Self;

  /// Superclass.
  typedef Optimizer Superclass;

  /** Constructor.
   * Hand functional and callback to parent class and initialize local 
   * variables.
   *@param stepFactor Factor by which the search step size is decreased.
   */
  BestNeighbourOptimizer ( const Self::ParameterType stepFactor = 0.5 )
  { 
    StepFactor = stepFactor; 
  };
 
  /** Perform the optimization.
   */
  virtual CallbackResult Optimize( CoordinateVector&, const Self::ParameterType = 1, const Self::ParameterType = 0 );

private:
  /** Search step factor.
   * This variable determines the factor by which to decrease the search step
   * size if no further improvement is possible at a certain resolution.
   * Reasonable values are in the range 0 < StepFactor < 1. For most cases,
   * a value of 0.5 has been found to provide optimum accuracy and performance.
   */
  Self::ParameterType StepFactor;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkBestNeighbourOptimizer_h_included_
