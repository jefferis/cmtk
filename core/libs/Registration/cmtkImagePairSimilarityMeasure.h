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

#ifndef __cmtkImagePairSimilarityMeasure_h_included_
#define __cmtkImagePairSimilarityMeasure_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkTypes.h>
#include <Base/cmtkFunctional.h>
#include <Base/cmtkUniformVolumeInterpolatorBase.h>

#include <System/cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base class for voxel metrics with pre-converted image data.
 */
class ImagePairSimilarityMeasure
{
public:
  /// This type.
  typedef ImagePairSimilarityMeasure Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Return type: same as cmtk::Functional.
  typedef Functional::ReturnType ReturnType;

  /** Constructor.
   */
  ImagePairSimilarityMeasure( const UniformVolume::SmartConstPtr& refVolume /*!< The reference image.*/,
			      const UniformVolume::SmartConstPtr& fltVolume /*!< The floating image.*/,
			      const Interpolators::InterpolationEnum interpolation = Interpolators::DEFAULT /*!< User-selected interpolation kernel*/ );

  /** Default constructor.
   */
  ImagePairSimilarityMeasure( const Interpolators::InterpolationEnum interpolation = Interpolators::DEFAULT /*!< User-selected interpolation kernel*/ ) : m_InterpolationMethod( interpolation ) {};

  /// Virtual destructor.
  virtual ~ImagePairSimilarityMeasure() {};

  /// Set reference volume.
  virtual void SetReferenceVolume( const UniformVolume::SmartConstPtr& refVolume );

  /** Set floating volume.
   * When the floating volume is set, a new interpolator object is also created.
   */
  virtual void SetFloatingVolume( const UniformVolume::SmartConstPtr& fltVolume );

  /// Reset metric computation.
  virtual void Reset() {}

  /// Get a value from the X distribution (reference image).
  Types::DataItem GetSampleX ( const size_t index ) const
  { 
    Types::DataItem data;
    this->m_ReferenceData->Get( data, index );
    return data;
  }

  /// Get a value from the X distribution (reference image).
  bool GetSampleX ( Types::DataItem& sample, const size_t index ) const
  { 
    return this->m_ReferenceData->Get( sample, index );
  }
  
  /// Get number of samples in the X data (reference image pixels).
  size_t GetNumberOfSamplesX() const
  {
    return this->m_ReferenceData->GetDataSize();
  }

  /// Get value range of X data (reference data).
  const Types::DataItemRange GetDataRangeX() const
  {
    return this->m_ReferenceData->GetRange();
  }
  
  /// Interpolate a value from the Y distribution (floating image).
  Types::DataItem GetSampleY( const int* index, const Types::Coordinate* frac ) const
  {
    return this->m_FloatingImageInterpolator->GetDataDirect( index, frac );
  }

  /// Get number of samples in the Y data (floating image pixels).
  size_t GetNumberOfSamplesY() const
  {
    return this->m_FloatingData->GetDataSize();
  }
  
  /// Get value range of Y data (floating data).
  const Types::DataItemRange GetDataRangeY() const
  {
    return this->m_FloatingData->GetRange();
  }
  
  /// Get scaled floating value if this metric rescales (implemented in derived classes), or input value if it does not (done here as the default).
  virtual Types::DataItem GetFloatingValueScaled( const Types::DataItem value ) const
  {
    return value;
  }

  /// Get the value of the metric.
  virtual Self::ReturnType Get() const = 0;

private:
  /// Smart pointer to reference volume.
  UniformVolume::SmartConstPtr m_ReferenceVolume;
  
  /// Smart pointer to reference image data.
  TypedArray::SmartConstPtr m_ReferenceData;
  
  /// Smart pointer to floating volume.
  UniformVolume::SmartConstPtr m_FloatingVolume;
  
  /// Smart pointer to floating image data.
  TypedArray::SmartConstPtr m_FloatingData;

  /// Interpolation method ID.
  Interpolators::InterpolationEnum m_InterpolationMethod;

  /// Floating image interpolator.
  cmtk::UniformVolumeInterpolatorBase::SmartConstPtr m_FloatingImageInterpolator;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSimilarityMeasure_h_included_
