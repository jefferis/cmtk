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

#ifndef __cmtkHistogramThresholdByVolume_h_included_
#define __cmtkHistogramThresholdByVolume_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkHistogram.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Class for computing binarization threshold from intensity histogram based on desired volume.
 */
template<class THistogram>
class HistogramThresholdByVolume
{
public:
  /// This class.
  typedef HistogramThresholdByVolume<THistogram> Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const for  this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// The histogram template parameter type.
  typedef THistogram HistogramType;

  /// Constructor: compute minimum threshold such that the number of samples above threshold is at least "volumeAbove".
  HistogramThresholdByVolume( const typename Self::HistogramType& histogram, const typename Self::HistogramType::BinType volumeAbove );

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

#include "cmtkHistogramThresholdByVolume.txx"

#endif // #ifndef __cmtkHistogramThresholdByVolume_h_included_
