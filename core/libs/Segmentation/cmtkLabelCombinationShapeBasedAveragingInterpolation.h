/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#ifndef __cmtkLabelCombinationShapeBasedAveragingInterpolation_h_included_
#define __cmtkLabelCombinationShapeBasedAveragingInterpolation_h_included_

#include <cmtkconfig.h>

#include <Segmentation/cmtkLabelCombinationShapeBasedAveraging.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkXformUniformVolume.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkSmartConstPtr.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Label image combination by Shape Based Averaging.
 *\see T. Rohlfing and C. R. Maurer, Jr., "Shape-based averaging," IEEE Transactions on Image Processing, vol. 16, no. 1, pp. 153-161, 2007. http://dx.doi.org/10.1109/TIP.2006.884936
 */
class LabelCombinationShapeBasedAveragingInterpolation 
  : private LabelCombinationShapeBasedAveraging
{
public:
  /// This class.
  typedef LabelCombinationShapeBasedAveragingInterpolation Self;

  /// Parent class.
  typedef LabelCombinationShapeBasedAveraging Superclass;
  
  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Label index type.
  typedef Superclass::LabelIndexType LabelIndexType;
  
  /// Real-value type for distance maps.
  typedef Superclass::DistanceMapRealType DistanceMapRealType;
  
  /// Constructor: compute label combination.
  LabelCombinationShapeBasedAveragingInterpolation( const std::vector<UniformVolume::SmartConstPtr>& labelImages /*!< Input label images. */,
						    const std::vector<cmtk::XformUniformVolume::SmartConstPtr> xformsToLabelImages /*!< Transformations with pre-assigned reference image grid.*/,
						    const UniformVolume::SmartConstPtr& targetGrid /*!< Target grid for all transformations. */,
						    const Self::LabelIndexType numberOfLabels = 0 /*!< Number of labels. If zero, the highest label index is determined from the data */ );
  
  /// Get result.
  TypedArray::SmartPtr GetResult() const;

private:
  /// Target grid for all transformations.
  const UniformVolume::SmartConstPtr m_TargetGrid;

  /// Vector of transformations to the label images.
  const std::vector<cmtk::XformUniformVolume::SmartConstPtr> m_Transformations;
};

} // namespace cmtk

#endif // #ifndef __cmtkLabelCombinationShapeBasedAveragingInterpolation_h_included_
