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

#ifndef __cmtkImagePairAffineRegistrationFunctionalTemplate_h_included_
#define __cmtkImagePairAffineRegistrationFunctionalTemplate_h_included_

#include <cmtkconfig.h>

#include <cmtkImagePairAffineRegistrationFunctional.h>

#include <cmtkInterpolator.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Functional that evaluates a voxel-based similarity measure.
 * This class defines the type of functional that is optimized during
 * voxel-based registration. It holds references to reference and floating data
 * and computes similarity as well as its gradient w.r.t. a given
 * transformation.
 *
 * The metric to be optimized is given by a template parameter, therefore 
 * allowing inlined code to be generated for efficient evaluation.
 */
template<class VM>
class ImagePairAffineRegistrationFunctionalTemplate :
  /// Inherit from affine voxel matching functional
  public ImagePairAffineRegistrationFunctional
{
public:
  /// This class type.
  typedef ImagePairAffineRegistrationFunctionalTemplate<VM> Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef ImagePairAffineRegistrationFunctional Superclass;

  /// Return type.
  typedef Functional::ReturnType ReturnType;

  /** Constructor.
   * Init pointers to volume and transformation objects and initialize
   * internal data structures.
   *@param reference The reference (i.e. static) volume.
   *@param floating The floating (i.e. transformed) volume.
   *@param xform A transformation template. This object determines the type
   * of transformation to be optimized. Its initial value is not relevant.
   */
  ImagePairAffineRegistrationFunctionalTemplate( UniformVolume::SmartPtr& reference, 
						 UniformVolume::SmartPtr& floating, 
						 const Interpolators::InterpolationEnum interpolation,
						 AffineXform::SmartPtr& affineXform )
    : ImagePairAffineRegistrationFunctional( reference, floating, affineXform )
  {
    this->m_NumberOfThreads = ThreadPool::GetGlobalThreadPool().GetNumberOfThreads();

    this->m_Metric = new VM( reference, floating, interpolation );
    this->m_ThreadMetric = Memory::AllocateArray<VM*>( m_NumberOfThreads );
    for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
      this->m_ThreadMetric[thread] = new VM( *(this->m_Metric) );
  }

  /// Destructor.
  virtual ~ImagePairAffineRegistrationFunctionalTemplate() 
  {
    for ( size_t thread = 0; thread < m_NumberOfThreads; ++thread )
      delete m_ThreadMetric[thread];
    Memory::DeleteArray( this->m_ThreadMetric );

    delete this->m_Metric;
  }

  /// Evaluate with new parameter vector.
  virtual typename Self::ReturnType EvaluateAt ( CoordinateVector& v ) 
  {
    this->m_AffineXform->SetParamVector( v );
    return this->Evaluate();
  }

  /** Compute functional value with volume clipping.
   * This function iterates over all voxels of the reference image that - after
   * applying the current coordinate transformation - are located inside the
   * mode image. This set of voxels is determined on-the-fly by an extension of
   * Liang and Barsky's "Parameterized Line-Clipping" technique.
   *
   * From the resulting sequence of reference/floating voxel pairs, the 
   * selected voxel-based similarity measure (metric) is computed.
   *@param v The current parameter vector describing the effective coordinate
   * transformation.
   *@return The computed similarity measure as returned by the "Metric" 
   * subobject.
   *@see VolumeClipping
   */
  virtual typename Self::ReturnType Evaluate();
  
  /** Number of threads that this object was created for.
   * This is the actual maximum number of threads running at any time, but not
   * necessarily the number of parallel tasks to be completed.
   * All duplicated data structures are generated with the multiplicity given
   * by this value. It is determined from Threads when the object is first
   * instanced. It cannot be changed afterwards.
   */
  size_t m_NumberOfThreads;

  /// Metric objects for the separate threads.
  VM** m_ThreadMetric;

  /// Mutex lock for access to global Metric field.
  MutexLock m_MetricMutex;

  /** Thread parameter block for incremental gradient computation.
   * This structure holds all thread-specific information. A pointer to an
   * instance of this structure is given to EvaluateGradientThread() for
   * each thread created.
   */
  typedef struct 
  {
    /// Pointer to the functional object that created the thread.
    Self *thisObject;
    /// Axes hash.
    const VolumeAxesHash* AxesHash;
    /// First plane of clipped reference volume.
    DataGrid::IndexType::ValueType StartZ;
    /// Last plane of clipped reference volume.
    DataGrid::IndexType::ValueType EndZ;
  } EvaluateTaskInfo;
 
  /// Info blocks for parallel threads evaluating functional gradient.
  std::vector<typename Self::EvaluateTaskInfo> m_EvaluateTaskInfo;

  /** Compute functional gradient as a thread.
    * This function (i.e., each thread) iterates over all parameters of the
    * current warp transformation. Among all active (i.e., not disabled)
    * parameters, it selects the ones that have an index with modulus
    * equal to the threads index when divided by the total number of threads.
    * For these parameters, the thread computes the partial derivative of the
    * functional by finite-difference approximation.
    */
  static void EvaluateThread( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t );

private:
  /// The metric (similarity measure) object.
  VM* m_Metric;
};

//@}

} // namespace cmtk

#include <cmtkImagePairAffineRegistrationFunctionalTemplate.txx>

#endif // #ifndef __cmtkImagePairAffineRegistrationFunctionalTemplate_h_included_
