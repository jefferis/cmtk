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

#ifndef __cmtkImagePairSimilarityMeasureNMI_h_included_
#define __cmtkImagePairSimilarityMeasureNMI_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkImagePairSimilarityJointHistogram.h>

#include <Base/cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

#ifdef _MSC_VER
#pragma warning (disable:4521)
#endif
/** Base class for voxel metrics with pre-converted image data.
 */
class ImagePairSimilarityMeasureNMI :
  /// Inherit generic image pair similarity class.
  public ImagePairSimilarityJointHistogram
{
public:
  /// This type.
  typedef ImagePairSimilarityMeasureNMI Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Return type: same as cmtk::Functional.
  typedef Functional::ReturnType ReturnType;

  /** Constructor.
   * For reference and model volume, InitDataset is called.
   *@param refVolume The reference (fixed) volume.
   *@param modVolume The model (transformed) volume.
   *@param initData If this flag is set (default), then the internal 
   * representation of the pixel data for both volumes is also created.
   * Derived classes may want to prevent this if they define their own
   * specific initialization (e.g., igsVoxelMatchingJointHistogram).
   */
  ImagePairSimilarityMeasureNMI( UniformVolume::SmartConstPtr& refVolume, UniformVolume::SmartConstPtr& fltVolume,
				 const Interpolators::InterpolationEnum interpolation = Interpolators::DEFAULT )
    : ImagePairSimilarityJointHistogram( refVolume, fltVolume, interpolation )
  {}

  /** Default constructor.
   */
  ImagePairSimilarityMeasureNMI() {};

  /// Get the value of the metric.
  virtual Self::ReturnType Get() const 
  {
    double HX, HY;
    
    this->m_JointHistogram.GetMarginalEntropies(HX,HY);
    const double HXY = this->m_JointHistogram.GetJointEntropy();
    
    return static_cast<Self::ReturnType>( (HX + HY) / HXY );
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSimilarityMeasureNMI_h_included_
