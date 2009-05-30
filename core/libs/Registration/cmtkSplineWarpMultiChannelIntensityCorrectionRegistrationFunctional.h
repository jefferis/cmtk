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

#ifndef __cmtkSplineWarpMultiChannelIntensityCorrectionRegistrationFunctional_h_included_
#define __cmtkSplineWarpMultiChannelIntensityCorrectionRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkMultiChannelRMIRegistrationFunctional.h>
#include <cmtkTemplateMultiChannelRegistrationFunctional.h>
#include <cmtkSplineWarpMultiChannelRegistrationFunctional.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Class for spline warp multi-channel registration functional. */
template<class TMetricFunctional = MultiChannelRMIRegistrationFunctional<> >
class SplineWarpMultiChannelIntensityCorrectionRegistrationFunctional :
  /** Inherit from multi-channel registration functional base class. */
  public SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional>
{
public:
  /** This class. */
  typedef SplineWarpMultiChannelIntensityCorrectionRegistrationFunctional<TMetricFunctional> Self;

  /** Smart pointer. */
  typedef SmartPointer<Self> SmartPtr;

  /** This class. */
  typedef SplineWarpMultiChannelRegistrationFunctional<TMetricFunctional> Superclass;

  /** Metric data class. */
  typedef typename TMetricFunctional::MetricData MetricData;

private:
  /** Continue metric computation and store reformatted floating channels for local recomputation. */
  virtual void ContinueMetricStoreReformatted( MetricData& metricData, const size_t rindex, const Vector3D& fvector );

  /** Continue metric computation and store reformatted floating channels for local recomputation. */
  virtual void ContinueMetric( MetricData& metricData, const size_t rindex, const Vector3D& fvector );
};

//@}

} // namespace cmtk

#include <cmtkSplineWarpMultiChannelIntensityCorrectionRegistrationFunctional.txx>

#endif // #ifndef __cmtkSplineWarpMultiChannelIntensityCorrectionRegistrationFunctional_h_included_
