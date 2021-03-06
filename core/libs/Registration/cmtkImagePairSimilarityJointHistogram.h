/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <Registration/cmtkImagePairSimilarityMeasure.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkFunctional.h>
#include <Base/cmtkJointHistogram.h>

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

  /// Parent class.
  typedef ImagePairSimilarityMeasure Superclass;

  /// Return type: same as cmtk::Functional.
  typedef Functional::ReturnType ReturnType;

  /** Constructor.
   * For reference and model volume, InitDataset is called.
   *\param refVolume The reference (fixed) volume.
   *\param fltVolume The floating (moving, transformed) volume.
   *\param interpolation ID of the interpolator to use for the floating image.
   */
  ImagePairSimilarityJointHistogram( UniformVolume::SmartConstPtr& refVolume, UniformVolume::SmartConstPtr& fltVolume, const Interpolators::InterpolationEnum interpolation = Interpolators::DEFAULT );

  /** Default constructor.
   */
  ImagePairSimilarityJointHistogram() {};

  /** Virtual destructor.
   */
  virtual ~ImagePairSimilarityJointHistogram() {};

  /** Set reference volume.
   * In addition to setting the reference volume via the base class, this function
   * also performs pre-scaling and parameter selection using Self::PrescaleData().
   * Afterwards the joint histogram size is re-allocated.
   */
  virtual void SetReferenceVolume( const UniformVolume::SmartConstPtr& refVolume );

  /** Set floating volume.
   * In addition to setting the floating volume via the base class, this function
   * also performs pre-scaling and parameter selection using Self::PrescaleData().
   * Afterwards the joint histogram size is re-allocated.
   */
  virtual void SetFloatingVolume( const UniformVolume::SmartConstPtr& fltVolume );

  /// Reset computation: clear joint histogram.
  virtual void Reset () 
  {
    this->m_JointHistogram.Reset();
  }

  /** Add a pair of values to the metric.
   */
  template<class T> void Increment( const T a, const T b )
  {
    this->m_JointHistogram.Increment( static_cast<size_t>( a ), std::max<size_t>( 0, std::min<size_t>( this->m_NumberOfBinsY-1, static_cast<size_t>( b ) ) ) );
  }

  /** Remove a pair of values from the metric.
   */
  template<class T> void Decrement( const T a, const T b )
  {
    this->m_JointHistogram.Decrement( static_cast<size_t>( a ), std::max<size_t>( 0, std::min<size_t>( this->m_NumberOfBinsY-1, static_cast<size_t>( b ) ) ) );
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

  /// Get scaled floating value if this metric rescales (implemented in derived classes), or input value if it does not (done here as the default).
  virtual Types::DataItem GetFloatingValueScaled( const Types::DataItem value ) const
  {
    return static_cast<Types::DataItem>( floor(this->m_ScaleFactorFloating*value+this->m_ScaleOffsetFloating) );
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
  UniformVolume::SmartPtr PrescaleData( const UniformVolume::SmartConstPtr& volume /*!< Input volume.*/,
					size_t* numberOfBins /*!< Output: number of bins that the histogram should allocate for the output volume.*/,
					Types::DataItem* scaleFactor /*!< Data scaling factor.*/,
					Types::DataItem* scaleOffset /*!< Data scaling offset.*/ );

  /// Store reference data rescaling offset.
  Types::DataItem m_ScaleOffsetReference;

  /// Store reference data rescaling factor.
  Types::DataItem m_ScaleFactorReference;

  /// Store floating data rescaling offset.
  Types::DataItem m_ScaleOffsetFloating;

  /// Store floating data rescaling factor.
  Types::DataItem m_ScaleFactorFloating;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSimilarityJointHistogram_h_included_
