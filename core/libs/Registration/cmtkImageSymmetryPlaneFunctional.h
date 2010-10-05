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

#ifndef __cmtkImageSymmetryPlaneFunctional_h_included_
#define __cmtkImageSymmetryPlaneFunctional_h_included_

#include <cmtkconfig.h>

#include "cmtkImageSymmetryPlaneFunctionalBase.h"

#include <Registration/cmtkImagePairSimilarityMeasure.h>
#include <Registration/cmtkImagePairSimilarityMeasureMSD.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Functional for finding a symmetry plane in 3-D volumes.
 */
class ImageSymmetryPlaneFunctional :
  /// Inherit functional interface.
  public ImageSymmetryPlaneFunctionalBase
{
public:
  /// This class.
  typedef ImageSymmetryPlaneFunctional Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass
  typedef ImageSymmetryPlaneFunctionalBase Superclass;

  /// Type of metric we're using.
  typedef ImagePairSimilarityMeasureMSD MetricType;
  
  /// Constructor.
  ImageSymmetryPlaneFunctional( UniformVolume::SmartConstPtr& volume );

  /// Constructor with value range limits.
  ImageSymmetryPlaneFunctional( UniformVolume::SmartConstPtr& volume, const Types::DataItemRange& valueRange );

  /// Destructor.
  virtual ~ImageSymmetryPlaneFunctional() {}

  /// Compute functional value.
  virtual Self::ReturnType Evaluate();

private:
  /// Image similarity measure.
  Self::MetricType::SmartPtr m_Metric;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageSymmetryPlaneFunctional_h_included_
