/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkSimpleLevelset_h_included_
#define __cmtkSimpleLevelset_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkUniformVolume.h"

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Class for computing a simple two-phase levelset evolution.
 */
class SimpleLevelset
{
public:
  /// Constructor.
  SimpleLevelset( UniformVolume::SmartConstPtr& volume ) : m_Volume( volume ) {}

  /// Set filter sigma parameter.
  void SetFilterSigma( const Types::Coordinate filterSigma )
  {
    this->m_FilterSigma = filterSigma;
  }

  /// Set evolution time delta.
  void SetTimeDelta( const Types::Coordinate timeDelta )
  {
    this->m_TimeDelta = timeDelta;
  }

  /// Set levelset threshold.
  void SetLevelsetThreshold( const Types::Coordinate levelsetThreshold )
  {
    this->m_LevelsetThreshold = levelsetThreshold;
  }
  
  /// Initialize levelset with a centered sphere.
  void InitializeCenteredSphere();

  /// Levelset evolution.
  virtual void Evolve( const int numberOfIterations /**!< Number of iterations */, 
		       const bool forceIterations = false /**!< If this is set, evolution continues until maximum iteration count is reached, even when convergence is detected */,
		       const bool verbose = false /**!< Verbose flag */ );

  /** Return levelset, optionally converting to a binarized byte pixel representation.
   *\warning If the levelset is retrieved with the "binarize" flag set, then the
   *  levelset stored in this object will remain binarized after the call and
   *  should be re-initialized before calling "Evolve" again.
   */
  UniformVolume::SmartPtr& GetLevelset( const bool binarize = false, /**!< If set, levelset is binarized and converted to byte data */ const float threshold = 0.5 /**!< Threshold for optional binarization */ );

protected:
  /// The volume to compute a levelset segmentation for.
  UniformVolume::SmartConstPtr m_Volume;

  /// The evolving levelset.
  UniformVolume::SmartPtr m_Levelset;

  /// Sigma parameter of the Gaussian filter kernel.
  Types::Coordinate m_FilterSigma;

  /// Delta time constant.
  Types::Coordinate m_TimeDelta;

  /// Levelset threshold.
  Types::Coordinate m_LevelsetThreshold;
};

} // namespace cmtk

#endif // #ifndef __cmtkSimpleLevelset_h_included_
