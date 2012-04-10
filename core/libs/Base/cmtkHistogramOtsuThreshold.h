/*
//
//  Copyright 2012 SRI International
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

#ifndef __cmtkHistogramOtsuThreshold_h_included_
#define __cmtkHistogramOtsuThreshold_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkHistogram.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Class for computing binarization threshold from intensity histogram using Otsu's method.
 *\see https://en.wikipedia.org/wiki/Otsu%27s_method
 *\see N. Otsu (1979). "A threshold selection method from gray-level histograms." IEEE Trans. Sys., Man., Cyber. 9 (1): 62–66. 
 * http://dx.doi.org/10.1109/TSMC.1979.4310076
 *\see D.-Y. Huang and C.-H. Wang, "Optimal multi-level thresholding using a two-stage Otsu optimization approach," Pattern Recognition Letters, vol. 30, no. 3, 
 * pp. 275-284, 2009. http://dx.doi.org/10.1016/j.patrec.2008.10.003
 */
template<class THistogram>
class HistogramOtsuThreshold
{
public:
  /// This class.
  typedef HistogramOtsuThreshold<THistogram> Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const for  this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// The histogram template parameter type.
  typedef THistogram HistogramType;

  /// Constructor: compute registrations.
  HistogramOtsuThreshold( const Self::HistogramType& histogram );

  /// Get the computed threshold.
  Types::DataItem Get() const 
  {
    return this->m_Threshold;
  }
  
private:
  /// Computed threshold.
  Types::DataItem m_Threshold;
};

//@}

} // namespace cmtk

#include "cmtkHistogramOtsuThreshold.txx"

#endif // #ifndef __cmtkHistogramOtsuThreshold_h_included_
