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

#ifndef __cmtkEntropyMinimizationIntensityCorrectionFunctionalBase_h_included_
#define __cmtkEntropyMinimizationIntensityCorrectionFunctionalBase_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkFunctional.h"
#include "System/cmtkSmartPtr.h"
#include "System/cmtkThreads.h"

#include "Base/cmtkUniformVolume.h"
#include "Base/cmtkDataGrid.h"
#include "Base/cmtkHistogram.h"
#include "Base/cmtkLogHistogram.h"
#include "Base/cmtkTemplateArray.h"

#include <vector>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{
/// Base class for entropy-minimzation MR bias correction functional.
class EntropyMinimizationIntensityCorrectionFunctionalBase
  : public Functional
{
public:
  /// This class type.
  typedef EntropyMinimizationIntensityCorrectionFunctionalBase Self;

  /// Pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass type.
  typedef Functional Superclass;

  /// Constructor.
  EntropyMinimizationIntensityCorrectionFunctionalBase() 
    : m_SamplingDensity( 1.0 ),
      m_NumberOfHistogramBins( 256 ),
      m_UseLogIntensities( false )
  {}

  /// Virtual destructor.
  virtual ~EntropyMinimizationIntensityCorrectionFunctionalBase() {}

  /// Get number of additive monomials.
  virtual size_t GetNumberOfMonomialsAdd() const = 0;

  /// Get number of multiplicative monomials.
  virtual size_t GetNumberOfMonomialsMul() const = 0; 

  /// Set input image.
  virtual void SetInputImage( UniformVolume::SmartConstPtr& inputImage );

  /// Set foreground mask.
  virtual void SetForegroundMask( const UniformVolume& foregroundMask );

  /// Set sampling density.
  virtual void SetSamplingDensity( const float samplingDensity )
  {
    this->m_SamplingDensity = samplingDensity;
  }

  /// Set number of histogram bins.
  virtual void SetNumberOfHistogramBins( const size_t numberOfHistogramBins )
  {
    this->m_NumberOfHistogramBins = numberOfHistogramBins;
  }

  /** Set flag for use of log intensities for entropy estimation.
   * Using log intensities compensates for the entropy increase otherwise caused by
   * spreading distributions of values in brightened areas. This can help make
   * bias field estimation more robust, potentially without any masking.
   */
  void SetUseLogIntensities( const bool flag )
  {
    this->m_UseLogIntensities = flag;
  }

  /// Get corrected output image.
  virtual UniformVolume::SmartPtr& GetOutputImage( const bool update = false )
  {
    if ( update )
      this->UpdateOutputImage( false /*foregroundOnly*/);
    
    return this->m_OutputImage;
  }

  /// Update and return corrected output image.
  virtual UniformVolume::SmartPtr& GetOutputImage( CoordinateVector& v, const bool foregroundOnly = false );
  
  /// Get additive bias field.
  virtual UniformVolume::SmartPtr GetBiasFieldAdd( const bool updateCompleteImage = false )
  {
    if ( updateCompleteImage )
      this->UpdateBiasFieldAdd( false /*foregroundOnly*/ );

    UniformVolume::SmartPtr biasField( this->m_OutputImage->CloneGrid() );
    biasField->SetData( this->m_BiasFieldAdd );
    return biasField;
  }

  /// Set additive bias field.
  virtual void SetBiasFieldAdd( const UniformVolume& biasFieldAdd )
  {
    biasFieldAdd.GetData()->BlockCopy( *(this->m_BiasFieldAdd), 0, 0, this->m_BiasFieldAdd->GetDataSize() );
  }

  /// Get multiplicative bias field.
  virtual UniformVolume::SmartPtr GetBiasFieldMul( const bool updateCompleteImage = false )
  {
    if ( updateCompleteImage )
      this->UpdateBiasFieldMul( false /*foregroundOnly*/ );

    UniformVolume::SmartPtr biasField( this->m_OutputImage->CloneGrid() );
    biasField->SetData( this->m_BiasFieldMul );
    return biasField;
  }

  /// Set multiplicative bias field.
  virtual void SetBiasFieldMul( const UniformVolume& biasFieldMul )
  {
    biasFieldMul.GetData()->BlockCopy( *(this->m_BiasFieldMul), 0, 0, this->m_BiasFieldMul->GetDataSize() );
  }

  /// Evaluate functional.
  virtual Self::ReturnType Evaluate()
  {
    return static_cast<Self::ReturnType>( -this->m_OutputImage->GetData()->GetEntropy( *this->m_EntropyHistogram ) );
  }

  /// Evaluate functional for given parameter vector.
  virtual Self::ReturnType EvaluateAt( CoordinateVector& v );

protected:
  /// Original input image.
  UniformVolume::SmartConstPtr m_InputImage;

  /// Input intensity image range.
  Types::DataItem m_InputImageRange;

  /// Evolving corrected output image.
  UniformVolume::SmartPtr m_OutputImage;

  /// Type for histogram.
  typedef Histogram<unsigned int> HistogramType;

  /// Type for histogram using log-intensities.
  typedef LogHistogram<unsigned int> LogHistogramType;

  /// Histogram for entropy evaluation.
  HistogramType::SmartPtr m_EntropyHistogram;

  /// Binary foreground mask.
  std::vector<bool> m_ForegroundMask;

  /// Update polynomial correction factors from input image.
  virtual void UpdateCorrectionFactors() = 0;

  /** Update output image estimate based on current bias field parameters.
   *\param allPixels If this flag is set, any existing image foreground mask is ignored
   * and all image pixels are updated. Default: only update foreground pixels.
   */
  virtual void UpdateOutputImage( const bool foregroundOnly = true );

  /// Additive bias field.
  FloatArray::SmartPtr m_BiasFieldAdd;

  /// Multiplicative bias field.
  FloatArray::SmartPtr m_BiasFieldMul;

  /// Jointly update both bias images.
  virtual void UpdateBiasFields( const bool foregroundOnly = true ) = 0;

  /// Update additive bias image.
  virtual void UpdateBiasFieldAdd( const bool foregroundOnly = true ) = 0;

  /// Update additive bias image.
  virtual void UpdateBiasFieldMul( const bool foregroundOnly = true ) = 0;

  /// Number of input image pixels.
  size_t m_NumberOfPixels;

protected:  
  /** Sampling density.
   * This defines the fraction of foreground pixels that are considered in
   * the computation.
   */
  float m_SamplingDensity;

  /// Number of histogram bins for entropy estimation.
  size_t m_NumberOfHistogramBins;

  /// Flag for using log-intensities for entropy estimation.
  bool m_UseLogIntensities;

private:  
  /// Class for output image update thread parameters.
  class UpdateOutputImageThreadParameters :
    public ThreadParameters<Self>
  {
  public:
    /// Flag as given to UpdateOutputImage().
    bool m_ForegroundOnly;
  };

  /// Thread function: update output image.
  static void UpdateOutputImageThreadFunc( void* args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkEntropyMinimizationIntensityCorrectionFunctionalBase_h_included_

