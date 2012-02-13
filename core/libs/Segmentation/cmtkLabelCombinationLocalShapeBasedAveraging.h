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

#ifndef __cmtkLabelCombinationLocalShapeBasedAveraging_h_included_
#define __cmtkLabelCombinationLocalShapeBasedAveraging_h_included_

#include <cmtkconfig.h>

#include <Segmentation/cmtkLabelCombinationLocalVoting.h>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Segmentation combination by locally-weighted shape-based averaging.
 */
class LabelCombinationLocalShapeBasedAveraging
  : public LabelCombinationLocalVoting
{
public:
  /// This class.
  typedef LabelCombinationLocalShapeBasedAveraging Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const to this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Parent class.
  typedef LabelCombinationLocalVoting Superclass;

  /// Constructor: compute label combination.
  LabelCombinationLocalShapeBasedAveraging( const UniformVolume::SmartConstPtr targetImage ) : Superclass( targetImage ), m_MaxLabelValue( 0 ) {}

  /// Set flag to detect local outliers at each pixel in the co-registered distance maps.
  void SetDetectLocalOutliers( const bool detectOutliers = true )
  {
    this->m_DetectLocalOutliers = detectOutliers;
  }
  
  /// Get resulting combined segmentation.
  virtual TypedArray::SmartPtr GetResult() const;  

private:
  /** Maximum label value.
   * If this is zero, then the class will determine the largest value that is present in the input label maps.
   */
  int m_MaxLabelValue;
  
  /// Flag for outlier detection.
  bool m_DetectLocalOutliers;

  /// Compute result for a region.
  void ComputeResultForRegion( TypedArray& result /*!< Evolving result label map */, 
			       std::vector<float>& resultDistance /*!< Evolving result distance map */, 
			       const int label /*!< Current label */,
			       const Self::TargetRegionType& region /*!< Image region to work on */, 
			       std::vector<UniformVolume::SmartConstPtr> dmaps /*!< Vector of distance maps per atlas for the current label */ ) const;
  
};

} // namespace cmtk

#endif // #ifndef __cmtkLabelCombinationLocalShapeBasedAveraging_h_included_
