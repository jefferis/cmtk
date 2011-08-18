/*
//
//  Copyright 2008-2011 SRI International
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

#include <Base/cmtkTypedArray.h>
#include <Base/cmtkHistogram.h>

namespace 
cmtk
{

/** Estimate noise level in data stored in a TypedArray.
 * Estimate Gaussian noise variance using naive peak finding method.
 *\author Torsten Rohlfing
 */
class
TypedArrayNoiseEstimatorNaiveGaussian
{
public:
  /// This class.
  typedef TypedArrayNoiseEstimatorNaiveGaussian Self;

  /// Constructor.
  TypedArrayNoiseEstimatorNaiveGaussian( const TypedArray& data, const size_t histogramBins = 255 );
  
  /// Get noise level.
  Types::DataItem GetNoiseLevelSigma() const
  {
    return this->m_NoiseLevelSigma;
  }

protected:
  /// Default constructor; should not be invoked by user code.
  TypedArrayNoiseEstimatorNaiveGaussian()
  {
    this->m_NoiseLevelSigma = 0;
  }

  /// The estimate noise sigma.
  Types::DataItem m_NoiseLevelSigma;
};

} // namespace cmtk
