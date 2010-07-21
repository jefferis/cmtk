/*
//
//  Copyright 2008-2009 SRI International
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

#include <cmtkconfig.h>

#include "Base/cmtkHistogram.h"
#include "Unstable/cmtkTypedArrayNoiseEstimatorBrummer.h"

namespace 
cmtk
{

/** Estimate noise level in data stored in a TypedArray.
 * Estimate Rician noise variance using Chang's maximum likelihood method.
 *\author Mike Hasak
 */
class
TypedArrayNoiseEstimatorMaximumLikelihood :
    /// Inherit basic fields and helper functions
    protected TypedArrayNoiseEstimatorBrummer
{
public:
  /// This class.
  typedef TypedArrayNoiseEstimatorMaximumLikelihood Self;

  /// Base class.
  typedef TypedArrayNoiseEstimatorBrummer Superclass;

  /// Constructor.
  TypedArrayNoiseEstimatorMaximumLikelihood( const TypedArray* data, const size_t histogramBins = 255 );
};

} // namespace cmtk
