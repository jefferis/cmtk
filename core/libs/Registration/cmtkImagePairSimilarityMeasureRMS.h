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

#ifndef __cmtkImagePairSimilarityMeasureRMS_h_included_
#define __cmtkImagePairSimilarityMeasureRMS_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkImagePairSimilarityMeasureMSD.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkTypedArray.h>
#include <System/cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{
#ifdef _MSC_VER
#pragma warning (disable:4521)
#endif
/** Mean squared difference metric.
 */
class ImagePairSimilarityMeasureRMS :
  /// Inherit MSD similarity measure.
  public ImagePairSimilarityMeasureMSD
{
public:
  /// This type.
  typedef ImagePairSimilarityMeasureRMS Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Parent class.
  typedef ImagePairSimilarityMeasureMSD Superclass;

  /** Constructor.
   * For reference and model volume, InitDataset is called.
   *\param refVolume The reference (fixed) volume.
   *\param fltVolume The floating (moving) volume.
   *\param interpolation ID of the interpolation algorithm to use for the floating image.
   */
  ImagePairSimilarityMeasureRMS( const UniformVolume::SmartConstPtr& refVolume, const UniformVolume::SmartConstPtr& fltVolume, const Interpolators::InterpolationEnum interpolation = Interpolators::DEFAULT );

  /** Virtual destructor.
   */
  virtual ~ImagePairSimilarityMeasureRMS() {};

  /// Get the value of the metric.
  virtual Self::ReturnType Get() const 
  {
    return -sqrt( -this->Superclass::Get() );
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSimilarityMeasureRMS_h_included_
