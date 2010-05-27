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

#ifndef __cmtkCongealingFunctional_h_included_
#define __cmtkCongealingFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkCongealingFunctionalBase.h>

#include <cmtkSmartPtr.h>
#include <cmtkThreads.h>

#include <cmtkUniformVolume.h>
#include <cmtkXform.h>
#include <cmtkHistogram.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Functional base class for groupwise congealing registration.
 * This functional evaluates Lilla Zollei's entropy criterion for massively
 * groupwise image registration.
 *
 *\References
 *
 * [1] L . Zoellei, E. Learned-Miller, E. Grimson, W.M. Wells III: "Efficient 
 *     Population Registration of 3D Data", ICCV 2005, Computer Vision for
 *     Biomedical Image Applications; Beijing, China
 */
template<class TXform>
class CongealingFunctional : 
  /** Inherit from template congealing base class. */
  public CongealingFunctionalBase<TXform>
{
public:
  /// Type of parent class.
  typedef CongealingFunctionalBase<TXform> Superclass;

  /// Type of this class.
  typedef CongealingFunctional<TXform> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Transformation type.
  typedef TXform XformType;

  /// Smart pointer to transformation type.
  typedef typename XformType::SmartPtr XformPointer;

  /// Base type for histogram bins.
  typedef unsigned int HistogramBinType;

  /// Histogram type.
  typedef Histogram<HistogramBinType> HistogramType;

  /// Constructor.
  CongealingFunctional();

  /// Destructor.
  virtual ~CongealingFunctional();

  /// Set number of histogram bins.
  virtual void SetNumberOfHistogramBins( const size_t numberOfHistogramBins );

  /** Set template grid.
   *\param The template grid that defines size and resolution for the
   *  implicit registration template.
   */
  virtual void SetTemplateGrid( UniformVolume::SmartPtr& templateGrid, const int downsample = 1, const bool useTemplateData = false );

  /// Evaluate functional with currently set parameters.
  virtual typename Self::ReturnType Evaluate();
  
protected:
  /// Standard deviation over all images by pixel.
  std::vector<byte> m_StandardDeviationByPixel;

#ifdef CMTK_BUILD_MPI
  /// Standard deviation over all images by pixel for current MPI node.
  std::vector<byte> m_StandardDeviationByPixelMPI;
#endif

  /// Update standard deviation by pixel.
  virtual void UpdateStandardDeviationByPixel();

  /// Flag whether standard deviations by pixel need updating.
  bool m_NeedsUpdateStandardDeviationByPixel;

  /** Create Gaussian kernels for samples in histogram.
   */
  void CreateGaussianKernels();

  /** Histogram sample kernels.
   * Each element is a pointer to an array that holds the elements of one side
   * of a discretely sampled symmetric kernel the represents samples with
   * uncertainty in a histogram. Element [0] of each such array is the central
   * element.
   */
  std::vector<HistogramBinType*> m_HistogramKernel;

  /** Radius of histogram sample kernel.
   * Each element here is the number of elements in the corresponding
   * HistogramKernel array.
   */
  std::vector<size_t> m_HistogramKernelRadius;

  /** Histograms for computation threads.
   */
  std::vector<HistogramType> m_ThreadHistograms;

  /// Update probabilistic sample table..
  virtual bool Wiggle();

private:
  /// Thread parameters with no further data.
  typedef ThreadParameters<Self> ThreadParametersType;

  /// Thread function to update standard dedviations by pixel.
  static void UpdateStandardDeviationByPixelThreadFunc( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t );
  
  /// Thread parameter for entropy evaluation.
  class EvaluateThreadParameters : 
    /// Inherit from generic thread parameter class.
    public ThreadParametersType
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
  static void EvaluateThread( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t threadCnt );

  /// Evaluate functional with currently set parameters with probabilistic sampling.
  static void EvaluateProbabilisticThread( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t threadCnt );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkCongealingFunctional_h_included_
