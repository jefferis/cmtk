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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkCongealingFunctionalBase_h_included_
#define __cmtkCongealingFunctionalBase_h_included_

#include <cmtkconfig.h>

#include <cmtkGroupwiseRegistrationFunctionalBase.h>

#include <cmtkSmartPtr.h>
#include <cmtkThreads.h>

#include <cmtkUniformVolume.h>
#include <cmtkXform.h>

#include <vector>

#ifdef CMTK_BUILD_MPI
#  include <mpi.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base class for groupwise registration functionals.
 */
template<class TXform,class THistogramBinType=unsigned int>
class CongealingFunctionalBase : 
  /** Inherit from generic groupwise functional. */
  public GroupwiseRegistrationFunctionalBase
{
public:
  /// Type of this class.
  typedef GroupwiseRegistrationFunctionalBase Superclass;
  
  /// Type of this class.
  typedef CongealingFunctionalBase<TXform,THistogramBinType> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Transformation type.
  typedef TXform XformType;

  /// Smart pointer to transformation type.
  typedef typename XformType::SmartPtr XformPointer;

  /// Base type for histogram bins.
  typedef THistogramBinType HistogramBinType;

  /// Constructor.
  CongealingFunctionalBase();

  /// Destructor.
  virtual ~CongealingFunctionalBase();

  /// Set number of histogram bins.
  virtual void SetNumberOfHistogramBins( const size_t numberOfHistogramBins );

  /// Set number of histogram bins.
  virtual void SetCropImageHistograms( const bool crop = true )
  {
    this->m_CropImageHistograms = crop;
  }

  /** Get coordinate transformation for one image in the group.
   *\param idx Index of the volume/transformation.
   *\return Transformation for the selected volume.
   */
  virtual const XformType* GetXformByIndex( const size_t idx ) const
  {
    return XformType::SmartPtr::DynamicCastFrom( this->m_XformVector[idx] );
  }

  /** Get coordinate transformation for one image in the group.
   *\param idx Index of the volume/transformation.
   *\return Transformation for the selected volume.
   */
  virtual typename XformType::SmartPtr GetXformByIndex( const size_t idx )
  {
    return XformType::SmartPtr::DynamicCastFrom( this->m_XformVector[idx] );
  }
  
  /** Get coordinate transformation for one active image in the group.
   *\param idx Index of the volume/transformation.
   *\return Transformation for the selected volume.
   */
  virtual const XformType* GetActiveXformByIndex( const size_t idx ) const
  {
    return XformType::SmartPtr::DynamicCastFrom( this->m_XformVector[idx + this->m_ActiveXformsFrom] );
  }

  /** Get coordinate transformation for one active image in the group.
   *\param idx Index of the volume/transformation.
   *\return Transformation for the selected volume.
   */
  virtual typename XformType::SmartPtr GetActiveXformByIndex( const size_t idx )
  {
    return XformType::SmartPtr::DynamicCastFrom( this->m_XformVector[idx + this->m_ActiveXformsFrom] );
  }

protected:
  /** Interpolate given moving image to template.
   *\param idx Index of of to reformat to template. This also determines which
   *  transformation is used.
   *\param destination The reformatted pixel data is stored in this array.
   *  Sufficient memory (for as many pixels as there are in the template grid)
   *  must be allocated there.
   */
  virtual void InterpolateImage( const size_t idx, byte* const destination );

private:
  /// Crop image histograms to get rid of high-intensity low-probability samples.
  bool m_CropImageHistograms;

  /// Thread parameters with no further data.
  typedef ThreadParameters<Self> ThreadParametersType;

  /// Thread function parameters for image interpolation.
  class InterpolateImageThreadParameters : 
    /// Inherit from generic thread parameters.
    public ThreadParametersType
  {
  public:
    /// Index of the image to be interpolated.
    size_t m_Idx;

    /// Pointer to storage that will hold the reformatted pixel data.
    byte* m_Destination;    
  };

  /// Image interpolation thread function.
  static CMTK_THREAD_RETURN_TYPE InterpolateImageThread( void* threadParameters );

  /// Image interpolation thread function with probabilistic sampling.
  static CMTK_THREAD_RETURN_TYPE InterpolateImageProbabilisticThread( void* threadParameters );

  /// Prepare data for one image.
  virtual UniformVolume* PrepareSingleImage( UniformVolume::SmartPtr& image );

  /// Smooth and pre-scale target images.
  virtual void PrepareTargetImages();

protected:
  /** User-defined background value from parent class, transformed to histogram bin index. */
  byte m_PrivateUserBackgroundValue;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkCongealingFunctionalBase_h_included_
