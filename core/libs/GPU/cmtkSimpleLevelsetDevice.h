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

#ifndef __cmtkSimpleLevelsetDevice_h_included_
#define __cmtkSimpleLevelsetDevice_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkUniformVolume.h"

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/** Class for computing a simple levelset evolution on the GPU
 */
class SimpleLevelsetDevice
{
public:
  /// Constructor.
  SimpleLevelsetDevice( UniformVolume::SmartConstPtr& volume ) : m_Volume( volume ) {}

  /// Set filter sigma parameter.
  void SetFilterSigma( const float filterSigma )
  {
    this->m_FilterSigma = filterSigma;
  }

  /// Set evolution time delta.
  void SetTimeDelta( const float timeDelta )
  {
    this->m_TimeDelta = timeDelta;
  }

  /// Set levelset threshold.
  void SetLevelsetThreshold( const float levelsetThreshold )
  {
    this->m_LevelsetThreshold = levelsetThreshold;
  }
  
  /// Initialize levelset with a centered sphere.
  void InitializeCenteredSphere();

  /// Levelset evolution.
  void Evolve( const int numberOfIterations /**!< Number of iterations */, const bool forceIterations = false /**!< If this is set, evolution continues until maximum iteration count is reached, even when convergence is detected */ );

  /// Return levelset, optionally converting to a binarized byte pixel representation.
  UniformVolume::SmartPtr& GetLevelset( const bool binarize = false, const float threshold = 0.5 );

private:
  /// The volume to compute a levelset segmentation for.
  UniformVolume::SmartConstPtr m_Volume;

  /// The evolving levelset.
  UniformVolume::SmartPtr m_Levelset;

  /// Sigma parameter of the Gaussian filter kernel.
  float m_FilterSigma;

  /// Delta time constant.
  float m_TimeDelta;

  /// Levelset threshold.
  float m_LevelsetThreshold;
};

} // namespace cmtk

#endif // #ifndef __cmtkSimpleLevelsetDevice_h_included_
