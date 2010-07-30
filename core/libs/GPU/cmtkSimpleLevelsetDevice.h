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

#include "Segmentation/cmtkSimpleLevelset.h"

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/** Class for computing a simple levelset evolution on the GPU
 */
class SimpleLevelsetDevice
    /// Inherit from CPU-based levelset class.
  : public SimpleLevelset
{
public:
  /// This class.
  typedef SimpleLevelsetDevice Self;

  /// Parent class.
  typedef SimpleLevelset Superclass;

  /// Constructor.
  SimpleLevelsetDevice( UniformVolume::SmartConstPtr& volume ) : Superclass( volume ) {}

  /// Levelset evolution on GPU.
  virtual void Evolve( const int numberOfIterations /**!< Number of iterations */, 
		       const bool forceIterations = false /**!< If this is set, evolution continues until maximum iteration count is reached, even when convergence is detected */ );
};

} // namespace cmtk

#endif // #ifndef __cmtkSimpleLevelsetDevice_h_included_
