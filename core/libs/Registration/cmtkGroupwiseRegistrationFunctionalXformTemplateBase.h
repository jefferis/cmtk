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

#ifndef __cmtkGroupwiseRegistrationFunctionalXformTemplateBase_h_included_
#define __cmtkGroupwiseRegistrationFunctionalXformTemplateBase_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkGroupwiseRegistrationFunctionalBase.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkThreads.h>
#include <System/cmtkThreadPool.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkXform.h>

#include <vector>

#ifdef CMTK_BUILD_MPI
#  include <mpi.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base class template for groupwise registration functionals.
 * This class template adds to its base class all basic functionality that depends
 * on the coordinate transformation model (affine vs. nonrigid) but does not require
 * implementation by explicit specialization. In other words, this class provides the
 * interface that is common to all transformation models.
 * 
 * The next level of derived classes exist in several specialized variants
 * that implement the transformation-dependent interfaces, i.e., member functions
 * that exist only for certain transformation models.
 */
template<class TXform>
class GroupwiseRegistrationFunctionalXformTemplateBase : 
  /** Inherit from generic groupwise functional. */
  public GroupwiseRegistrationFunctionalBase
{
public:
  /// Type of this class.
  typedef GroupwiseRegistrationFunctionalBase Superclass;
  
  /// Type of this class.
  typedef GroupwiseRegistrationFunctionalXformTemplateBase<TXform> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Transformation type.
  typedef TXform XformType;

  /// Smart pointer to transformation type.
  typedef typename XformType::SmartPtr XformPointer;

  /// Constructor.
  GroupwiseRegistrationFunctionalXformTemplateBase();

  /// Destructor.
  virtual ~GroupwiseRegistrationFunctionalXformTemplateBase();

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
  /** Number of (usable) histogram bins.
   */
  size_t m_HistogramBins;

  /** Maximal radius of histogram kernels.
   */
  size_t m_HistogramKernelRadiusMax;

  /** Threshold for maximum fraction of reformatted pixels from any given image that may be outside FOV.
   * If the number of outside pixels for any one image exceeds this threshold (as a fraction of
   * total number of reformatted pixels) then an exception is thrown.
   */
  float m_MaxRelativeNumberOutsidePixels;

  /** User-defined background value from parent class, transformed to histogram bin index. */
  byte m_PrivateUserBackgroundValue;

private:
  /// Crop image histograms to get rid of high-intensity low-probability samples.
  bool m_CropImageHistograms;

  /// Prepare data for one image.
  virtual UniformVolume::SmartPtr PrepareSingleImage( UniformVolume::SmartPtr& image );

  /// Smooth and pre-scale target images.
  virtual void PrepareTargetImages();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkGroupwiseRegistrationFunctionalXformTemplateBase_h_included_
