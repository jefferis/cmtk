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

#ifndef __cmtkSplineWarpCongealingFunctional_h_included_
#define __cmtkSplineWarpCongealingFunctional_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkCongealingFunctional.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkThreads.h>
#include <System/cmtkMutexLock.h>
#include <System/cmtkThreadSemaphore.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkHistogram.h>

#include <IO/cmtkClassStream.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Functional for spline warp congealing.
 * This functional evaluates Lilla Zollei's entropy criterion for massively
 * groupwise image registration.
 *
 *\section refs References
 *
 * [1] L . Zoellei, E. Learned-Miller, E. Grimson, W.M. Wells III: "Efficient 
 *     Population Registration of 3D Data", ICCV 2005, Computer Vision for
 *     Biomedical Image Applications; Beijing, China
 */
class SplineWarpCongealingFunctional : 
  public CongealingFunctional<SplineWarpXform>
{
public:
  /// Type of parent class.
  typedef CongealingFunctional<SplineWarpXform> Superclass;

  /// Superclass histogram type.
  typedef Superclass::HistogramType HistogramType;

  /// Type of this class.
  typedef SplineWarpCongealingFunctional Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /** Initialize spline warp transformations.
   */
  virtual void InitializeXformsFromAffine( const Types::Coordinate gridSpacing /*!< Control point grid spacing in real-world units*/,
					   std::vector<AffineXform::SmartPtr> initialAffineXformsVector /*!< Vector of initial affine coordinate transformations*/,
					   const bool exactSpacing = true /*!< If set, the control point spacing will be exactly as given in the first parameter*/ );

  /// Refine transformation control point grids.
  virtual void RefineTransformationGrids();

  /// Call inherited function and allocate local storage.
  virtual void SetTemplateGrid( UniformVolume::SmartPtr& templateGrid, const int downsample = 1, const bool useTemplateData = false );
    
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
  
  /// Update standard deviation by pixel.
  virtual void UpdateStandardDeviationByPixel();

  /// Entropies over all images by pixel for fast local recomputation.
  std::vector<double> m_EntropyByPixel;

#ifdef CMTK_USE_MPI
  /// Temporary storage for values computed on current MPI node.
  std::vector<double> m_EntropyByPixelMPI;
#endif

  /// Thread parameter for entropy evaluation.
  class EvaluateThreadParameters : 
    /// Inherit from generic thread parameter class.
    public ThreadParameters<Self>
  {
  public:
    /// Upon return from the thread function, this holds the partial entropy.
    double m_Entropy;

    /** Upon return from the thread function, this holds the number of
      * pixels with full image count, i.e., pixels that are within all
      * target images.
      */
    unsigned int m_Count;
  };
  
  /// Evaluate functional with currently set parameters.
  static void EvaluateThread( void* args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

  /// Thread function parameters for image interpolation.
  class EvaluateLocalGradientThreadParameters : 
    /// Inherit from generic thread parameters.
    public ThreadParameters<Self>
  {
  public:
    /// Global parameter step scale factor.
    Types::Coordinate m_Step;

    /// Pointer to array with computed local gradient components.
    Types::Coordinate* m_Gradient;

#ifdef CMTK_USE_MPI
    /// Index of first control point to be computed by all threads in this iteration.
    size_t m_FirstIndexToCompute;

#endif
  };

  /// Index of next control point to be processed by an available thread.
  size_t m_ControlPointIndexNext;

  /// Index of (last control point + 1) to be processed by an available thread.
  size_t m_ControlPointIndexLast;
  
  /// Mutex lock for control point queue.
  MutexLock m_ControlPointIndexLock;

  /// Class for static thread storage to avoid recurring memory allocations.
  class StaticThreadStorage
  {
  public:
    /// Initialize thread storage based on parent functional ("This").
    void Initialize( const Self* This );
    
    /// Function values evaluated at x+delta.
    std::vector<Self::ReturnType> m_FPlus;

    /// Function values evaluated at x-delta.
    std::vector<Self::ReturnType> m_FMinus;

    /// Pixel count for transformation +delta.
    std::vector<unsigned int> m_CountByParameterPlus;

    /// Pixel count for transformation -delta.
    std::vector<unsigned int> m_CountByParameterMinus;

    /// Copies of transformation objects.
    std::vector<SplineWarpXform::SmartPtr> m_Xforms;

    /// List of transformed vectors.
    std::vector<Vector3D> m_VectorList;
    
    /// Floating pixel count per template pixel.
    std::vector<size_t> m_Count;
    
    /// Stack histograms per pixel.
    std::vector<HistogramType> m_Histogram;

    /// Flag: do transformation parameters need to be copied from functional?
    bool m_NeedToCopyXformParameters;
  };

  /// Static thread storage array.
  std::vector<StaticThreadStorage> m_StaticThreadStorage;

#ifdef CMTK_USE_MPI
  /** Task function: Compute local gradient of the cost function for gradient approximation.
   * This function takes into consideration that in a spline warp, each control point
   * effects only a local neighborhood. It also groups the parameters by control
   * point and works over all images and x,y,z to speed things up substantially.
   */
  static CMTK_THREAD_RETURN_TYPE EvaluateLocalGradientThreadFunc( void* args );

  /// Enumeration type for MPI message tags.
  typedef enum {
    MESSAGE_TAG_COMPUTE = 1,
    MESSAGE_TAG_RESULTS = 2,
    MESSAGE_TAG_FINISHED = 3
  } MessageTagMPI;

  /// Thread work queue condition variable.
  ThreadSemaphore m_ThreadWorkSemaphore;

  /// Thread available condition variable.
  ThreadSemaphore m_ThreadReadySemaphore;

  /// Reorder gradient components received from other nodes into final gradient vector.
  void ReorderGradientComponents( Types::Coordinate *const dst, const Types::Coordinate* src, const size_t fromCpIdx, const size_t toCpIdx );
#else // #ifdef CMTK_USE_MPI
  /** Task function: Compute local gradient of the cost function for gradient approximation.
   * This function takes into consideration that in a spline warp, each control point
   * effects only a local neighborhood. It also groups the parameters by control
   * point and works over all images and x,y,z to speed things up substantially.
   */
  static void EvaluateLocalGradientThreadFunc( void* args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );
#endif
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSplineWarpCongealingFunctional_h_included_
