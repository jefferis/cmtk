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

#ifndef __cmtkImagePairSimilarityMeasureMI_h_included_
#define __cmtkImagePairSimilarityMeasureMI_h_included_

#include <cmtkconfig.h>

#include <cmtkImagePairSimilarityJointHistogram.h>

#include <cmtkUniformVolume.h>

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
class ImagePairSimilarityMeasureMI :
  /// Inherit generic image pair similarity class.
  public ImagePairSimilarityJointHistogram
{
public:
  /// This type.
  typedef ImagePairSimilarityMeasureMI Self;

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
  ImagePairSimilarityMeasureMI( const UniformVolume::SmartPtr& refVolume, const UniformVolume::SmartPtr& fltVolume,
				const Interpolators::InterpolationEnum interpolation = Interpolators::DEFAULT )
    : ImagePairSimilarityJointHistogram( refVolume, fltVolume, interpolation )
  {}

  /** Default constructor.
   */
  ImagePairSimilarityMeasureMI() {};

  /// Copy constructor.
  ImagePairSimilarityMeasureMI( const Self& src )
    : ImagePairSimilarityJointHistogram( src )
  {}

  /// Copy constructor.
  ImagePairSimilarityMeasureMI( Self& other, const bool copyData = false )
  : ImagePairSimilarityJointHistogram( other, copyData )
  {}

  /// Get the value of the metric.
  Self::ReturnType Get() const 
  {
    double HX, HY, HXY;
    
    this->m_JointHistogram.GetMarginalEntropies(HX,HY);
    this->m_JointHistogram.GetJointEntropy(HXY);
    
    return static_cast<Self::ReturnType>( HX + HY - HXY );
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSimilarityMeasureMI_h_included_
