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

#ifndef __cmtkBestDirectionOptimizer_h_included_
#define __cmtkBestDirectionOptimizer_h_included_

#include <cmtkconfig.h>

#include <cmtkOptimizer.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Optimizer derived from BestNeighbourOptimizer.
 */
class BestDirectionOptimizer : 
  /// Inherit generic optimizer features.
  public Optimizer 
{
public:
  /// This class.
  typedef BestDirectionOptimizer Self;

  /// Superclass.
  typedef Optimizer Superclass;

  /// Flag whether to use maximum (1) or Euclid (0) for normalization.
  cmtkGetSetMacroDefault(bool,UseMaxNorm,true);

  /** Threshold for direction components.
   * Before searching in a certain directions, all components below this
   * fraction of the maximum absolute component are set to zero. This is
   * done to prevent changes in irrelevant parameters and thus support local
   * recomputation when optimizing elastic deformation parameters for example.
   * set this flag to 0 to disable thresholding; set it to 1 to remove all but
   * the most significant components.
   */
  cmtkGetSetMacroDefault(Self::ParameterType,DirectionThreshold,-1);

  /** Number of repetitions of each search level, even if previously unsuccessful.
   * This is one by default, so whenever no successful update was made, the level
   * is finished. Setting this to values larger than 1 only makes sense if the
   * optimized functional is changing over time, e.g., due to probabilistic 
   * effects. 
   */
  cmtkGetSetMacro(int,RepeatLevelCount);

  /** Agressive mode.
   * If this flag is set, the optimization is continued at one level as long as
   * there is an improvement to the target function at any step level. Otherwise,
   * steps in the binary seach phase are not considered and search terminates 
   * earlier.
   */
  cmtkGetSetMacro(bool,AggressiveMode);

  /// Constructor.
  BestDirectionOptimizer ( const Self::ParameterType stepFactor = 0.5, const Self::ParameterType = 0.1 )
  { 
    StepFactor = stepFactor;
    this->m_UseMaxNorm = true;
    this->m_DirectionThreshold = -1;
    this->m_RepeatLevelCount = 1;
    this->m_AggressiveMode = false;
  };
 
  /// Optimize functional.
  virtual CallbackResult Optimize( CoordinateVector&, const Self::ParameterType = 1, const Self::ParameterType = 0 );

private:
  /// Factor by which the step size is reduced after each pass.
  Self::ParameterType StepFactor;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkBestDirectionOptimizer_h_included_
