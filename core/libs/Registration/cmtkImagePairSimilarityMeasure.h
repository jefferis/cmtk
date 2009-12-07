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

#ifndef __cmtkImagePairSimilarityMeasure_h_included_
#define __cmtkImagePairSimilarityMeasure_h_included_

#include <cmtkconfig.h>

#include <cmtkUniformVolume.h>
#include <cmtkTypes.h>
#include <cmtkSmartPtr.h>

#include <cmtkFunctional.h>

#include <cmtkUniformVolumeInterpolatorBase.h>

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
  ImagePairSimilarityMeasure
  ( const UniformVolume::SmartPtr& refVolume, //!< The reference image.
    const UniformVolume::SmartPtr& fltVolume, //!< The floating image.
    const Interpolators::InterpolationEnum interpolation = Interpolators::DEFAULT //!< User-selected interpolation kernel
    );

  /** Default constructor.
   */
  ImagePairSimilarityMeasure() {};

  /// Copy constructor.
  ImagePairSimilarityMeasure( const Self& src );

  /// Copy constructor.
  ImagePairSimilarityMeasure( Self& other, const bool copyData = false );

  /// Get a value from the X distribution (reference image).
  Types::DataItem GetSampleX ( const size_t index ) const
  { 
    Types::DataItem data;
    this->m_ReferenceData->Get( data, index );
    return data;
  }
  
  /// Get a value from the Y distribution (floating image).
  Types::DataItem GetSampleY ( const size_t index ) const 
  { 
    Types::DataItem data;
    this->m_FloatingData->Get( data, index );
    return data;
  }
  
  /// Interpolate a value from the Y distribution (floating image).
  Types::DataItem GetSampleY ( const size_t baseIndex, const int* index, const Types::Coordinate* frac ) const
  {
    return this->m_FloatingImageInterpolator->GetDataDirect( baseIndex, index, frac );
  }

private:
  /// Smart pointer to reference image data.
  const TypedArray::SmartPtr m_ReferenceData;
  
  /// Smart pointer to floating image data.
  const TypedArray::SmartPtr m_FloatingData;

  /// Floating image interpolator.
  const cmtk::UniformVolumeInterpolatorBase::SmartPtr m_FloatingImageInterpolator;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSimilarityMeasure_h_included_
