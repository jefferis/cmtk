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

#ifndef __cmtkImagePairSimilarityJointHistogram_h_included_
#define __cmtkImagePairSimilarityJointHistogram_h_included_

#include <cmtkconfig.h>

#include <cmtkImagePairSimilarityMeasure.h>

#include <cmtkUniformVolume.h>
#include <cmtkFunctional.h>
#include <cmtkJointHistogram.h>

#include <algorithm>

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
class ImagePairSimilarityJointHistogram :
  /// Inherit generic image pair similarity class.
  public ImagePairSimilarityMeasure
{
public:
  /// This type.
  typedef ImagePairSimilarityJointHistogram Self;

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
  ImagePairSimilarityJointHistogram( const UniformVolume::SmartPtr& refVolume, const UniformVolume::SmartPtr& fltVolume,
				     const Interpolators::InterpolationEnum interpolation = Interpolators::DEFAULT );

  /** Default constructor.
   */
  ImagePairSimilarityJointHistogram() {};

  /// Copy constructor.
  ImagePairSimilarityJointHistogram( const Self& src );

  /// Copy constructor.
  ImagePairSimilarityJointHistogram( Self& other, const bool copyData = false );

  /// Reset computation: clear joint histogram.
  void Reset () 
  {
    this->m_JointHistogram.Reset();
  }

  /** Add a pair of values to the metric.
   */
  template<class T> void Increment( const T a, const T b )
  {
    this->m_JointHistogram.Increment( static_cast<size_t>( a ), std::max<size_t>( 0, std::min<size_t>( this->m_NumberOfBinsY-1, static_cast<size_t>( b ) ) ) );
  }

  /// Add another metric object to this one.
  void Add ( const Self& other )
  {
    this->m_JointHistogram.AddJointHistogram( other.m_JointHistogram );
  }

  /// Add another metric object to this one.
  void Remove ( const Self& other )
  {
    this->m_JointHistogram.RemoveJointHistogram( other.m_JointHistogram );
  }

protected:
  /// Number of X bins (reference image)
  size_t m_NumberOfBinsX;

  /// Number of Y bins (floating image)
  size_t m_NumberOfBinsY;

  /// The joint histogram.
  JointHistogram<unsigned int> m_JointHistogram;

private:
  /** Duplicate and pre-scale image data so that we have the histogram bin numbers readily available.
   *\return A new volume with the same geometry as the input volume, but for DATACLASS_GREY, all pixel
   * values will have been rescaled to represent histogram bin indexes directly.
   */
  UniformVolume::SmartPtr PrescaleData
  ( const UniformVolume::SmartPtr& volume, //!< Input volume.
    size_t* numberOfBins //!< Output: number of bins that the histogram should allocate for the output volume.
    );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSimilarityJointHistogram_h_included_
