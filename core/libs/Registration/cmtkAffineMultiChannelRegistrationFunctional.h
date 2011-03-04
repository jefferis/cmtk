/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#ifndef __cmtkAffineMultiChannelRegistrationFunctional_h_included_
#define __cmtkAffineMultiChannelRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkTemplateMultiChannelRegistrationFunctional.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkVolumeClipping.h>
#include <Base/cmtkTransformedVolumeAxes.h>

#include <System/cmtkThreads.h>
#include <System/cmtkMutexLock.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Forward declaration of class template. */
template<class MetricFunctionalClass> 
class SplineWarpMultiChannelRegistrationFunctional;

/** Class for affine multi-channel registration functional. */
template<class TMultiChannelMetricFunctional>
class AffineMultiChannelRegistrationFunctional :
  /** Inherit from multi-channel registration functional base class. */
  public TemplateMultiChannelRegistrationFunctional<AffineXform,TMultiChannelMetricFunctional>
{
public:
  /** This class. */
  typedef AffineMultiChannelRegistrationFunctional<TMultiChannelMetricFunctional> Self;

  /** Smart pointer. */
  typedef SmartPointer<Self> SmartPtr;

  /** This class. */
  typedef TemplateMultiChannelRegistrationFunctional<AffineXform,TMultiChannelMetricFunctional> Superclass;
  
  /** Default constructor. */
  AffineMultiChannelRegistrationFunctional()
    : m_NumberOfThreads( Threads::GetNumberOfThreads() )
  {
  }

  /** Initialize transformation. */
  void InitTransformation( const bool alignCenters );

  /** Compute functional value with volume clipping using multi-threading.
   * This function iterates over all voxels of the reference image that - after
   * applying the current coordinate transformation - are located inside the
   * mode image. This set of voxels is determined on-the-fly by an extension of
   * Liang and Barsky's "Parameterized Line-Clipping" technique.
   *
   * From the resulting sequence of reference/floating voxel pairs, the 
   * selected voxel-based similarity measure (metric) is computed.
   *\return The computed similarity measure as returned by the "Metric" 
   * subobject.
   *\see VolumeClipping
   */
  virtual typename Self::ReturnType Evaluate();

private:
  /** Volume clipping object. */
  VolumeClipping m_VolumeClipper;

  /** Perform clipping/cropping in z-direction.
   * This function computes the intersection of reference and floating data in
   * z-direction. It determines the range of indices of those planes in the
   * reference that intersect the floating. This is the range over which to 
   * for-loop during metric computation.
   *\param clipper A volume clipping object with clipping boundaries and grid
   * orientation set.
   *\param origin Starting point of the reference volume.
   *\param start Upon return, this reference is set to the index of first plane
   * in the reference that intersects the floating.
   *\param end Upon return, this reference is set to one plus the index of the
   * last plane in the reference that intersects the floating.
   *\return true if there is an intersection of reference and floating, false if there
   * isn't. The range of indices returned in "start" and "end" is only
   * guaranteed to be valid if 1 is the return value.
   */
  bool ClipZ ( const VolumeClipping& clipper, const Vector3D& origin, int& start, int& end ) const;
    
  /** Perform clipping/cropping in x-direction.
   * This function computes the intersection of reference and floating data in
   * x-direction. It determines the range of indices of those voxels in the
   * current reference row that intersect the floating image. This is the range
   * over which to for-loop during metric computation.
   *
   * Compared to ClipZ and ClipY, this step has to operate very exact as there
   * is no further level that would reduce remaining invalid voxels. Therefore,
   * clipper.ClipX() is called with an extended initial range of indices and an
   * explicitly open upper bound.
   *
   * This is necessary to discriminate inside-boundary from on-boundary voxels.
   * For the right, upper and back boundary, on-boundary voxels are already
   * outside the allowed range as the upper boundaries of the volume are open
   * in terms of interpolation.
   *\param clipper A volume clipping object with clipping boundaries and grid
   * orientation set.
   *\param origin Starting point of the current row in the reference volume.
   *\param start Upon return, this reference is set to the index of first voxel
   * in the reference that intersects the floating image.
   *\param end Upon return, this reference is set to one plus the index of the
   * last voxel in the reference that intersects the floating image.
   *\return true if there is an intersection of the current reference row and
   * the floating, false if there isn't. The range of indices returned in "start"
   * and "end" is only guaranteed to be valid if true is the return value.
   */
  bool ClipX( const VolumeClipping& clipper, const Vector3D& origin, int& start, int &end ) const;

  /** Perform clipping/cropping in y-direction.
   * This function computes the intersection of reference and floating data in
   * y-direction. It determines the range of indices of those rows in the
   * current reference plane that intersect the floating image. This is the
   * range over which to for-loop during metric computation.
   *\param clipper A volume clipping object with clipping boundaries and grid
   * orientation set.
   *\param origin Starting point of the current plane in the reference volume.
   *\param start Upon return, this reference is set to the index of first row
   * in the reference that intersects the floating image.
   *\param end Upon return, this reference is set to one plus the index of the
   * last row in the reference that intersects the floating image.
   *\return true if there is an intersection of the current reference plane and
   * the floating, false if there isn't. The range of indices returned in "start" 
   * and "end" is only guaranteed to be valid if true is the return value.
   */
  bool ClipY( const VolumeClipping& clipper, const Vector3D& origin, int& start, int& end ) const;

  /** Number of parallel threads. */
  const size_t m_NumberOfThreads;

  /** Parameters for threaded metric computation. */
  class EvaluateThreadParameters : 
    /// Inherit from generic thread parameters.
    public ThreadParameters<Self>
  {
  public:
    /** Pointer to transformed reference axes arrays. */
    const TransformedVolumeAxes* m_TransformedAxes;

    /** First clipped pixel index in z direction. */
    int m_StartZ;

    /** Last clipped pixel index in z direction plus one. */
    int m_EndZ;
  };

  /** Thread function for metric computation. */
  static void EvaluateThreadFunction( void* args, const size_t taskIdx, const size_t taskCnt );

#ifdef CMTK_BUILD_SMP
  /** Mutex lock for shared metric data object. */
  MutexLock m_MetricDataMutex;
#endif

  /** Make spline functional friend to access this class' channel image vectors. */
  template<class MetricFunctionalClass> friend class SplineWarpMultiChannelRegistrationFunctional;
};

//@}

} // namespace cmtk

#include "cmtkAffineMultiChannelRegistrationFunctional.txx"

#endif // #ifndef __cmtkAffineMultiChannelRegistrationFunctional_h_included_
