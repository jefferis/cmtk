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

#ifndef __cmtkSplineWarpGroupwiseRegistrationRMIFunctional_h_included_
#define __cmtkSplineWarpGroupwiseRegistrationRMIFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkGroupwiseRegistrationRMIFunctional.h>

#include <cmtkSmartPtr.h>
#include <cmtkThreads.h>

#include <cmtkUniformVolume.h>
#include <cmtkSplineWarpXform.h>
#include <cmtkHistogram.h>

#include <vector>

#include <cmtkClassStream.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Functional for spline warp groupwise registration.
 */
class SplineWarpGroupwiseRegistrationRMIFunctional : 
  public GroupwiseRegistrationRMIFunctional<SplineWarpXform>
{
public:
  /// Type of parent class.
  typedef GroupwiseRegistrationRMIFunctional<SplineWarpXform> Superclass;

  /// Type of this class.
  typedef SplineWarpGroupwiseRegistrationRMIFunctional Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Evaluate functional with currently set parameters.
  virtual Self::ReturnType Evaluate();

  /// Evaluate functional and set parameters.
  virtual Self::ReturnType EvaluateAt( CoordinateVector& v )
  {
    return this->Superclass::EvaluateAt( v );
  }

  /** Compute functional value and gradient.
   *\param v Parameter vector.
   *\param g The extimated gradient of the functional is stored in this vector.
   *\param step Step size for finite difference gradient approximation. Default
   *  is 1 mm.
   *\return Const function value for given parameters.
   */
  virtual Self::ReturnType EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step = 1 );

private:
  /// Update deactivated control points.
  virtual void UpdateActiveControlPoints();
  
  /// Local information measure for neighborhood of each control point.
  std::vector<byte> m_InformationByControlPoint;

  /// Flag whether information by control point needs to be updated.
  bool m_NeedsUpdateInformationByControlPoint;

  /// Update local information by control point.
  virtual void UpdateInformationByControlPoint();

  /// Update control point schedule for gradient approximation.
  virtual void UpdateControlPointSchedule();

  /// Processing schedule for overlap-free parallel processing of control points.
  std::vector<int> m_ControlPointSchedule;

  /// Maximum number of concurrent jobs working on warps that is still overlap-free.
  size_t m_ControlPointScheduleOverlapFreeMaxLength;

  /// Thread function parameters for image interpolation.
  class EvaluateLocalGradientThreadParameters : 
    /// Inherit from generic thread parameters.
    public ThreadParameters<Self>
  {
  public:
    /// Unique thread storage index.
    size_t m_ThreadStorageIndex;

    /// Global parameter step scale factor.
    Types::Coordinate m_Step;

    /// Pointer to array with computed local gradient components.
    Types::Coordinate* m_Gradient;
    
    /// Current metric value.
    Self::ReturnType m_MetricBaseValue;

#ifdef CMTK_BUILD_MPI
    /// Index of first control point to be computed by all threads in this iteration.
    size_t m_FirstIndexToCompute;
#endif
  };

  /** Thread function: Compute local gradient of the cost function for gradient approximation.
   * This function takes into consideration that in a spline warp, each control point
   * effects only a local neighborhood. It also groups the parameters by control
   * point and works over all images and x,y,z to speed things up substantially.
   */
  static CMTK_THREAD_RETURN_TYPE EvaluateLocalGradientThreadFunc( void* args );

#ifdef CMTK_BUILD_MPI
  void ReorderGradientComponents( Types::Coordinate *const dst, const Types::Coordinate* src, const size_t fromCpIdx, const size_t toCpIdx );
#endif
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSplineWarpGroupwiseRegistrationRMIFunctional_h_included_
