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

#ifndef __cmtkImagePairSimilarityMeasureMSD_h_included_
#define __cmtkImagePairSimilarityMeasureMSD_h_included_

#include <cmtkconfig.h>

#include "Registration/cmtkImagePairSimilarityMeasure.h"

#include "Base/cmtkUniformVolume.h"
#include "Base/cmtkTypedArray.h"
#include "Base/cmtkMathUtil.h"
#include "System/cmtkSmartPtr.h"

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
class ImagePairSimilarityMeasureMSD :
  /// Inherit generic pairwise similarity measure class
  public ImagePairSimilarityMeasure
{
public:
  /// This type.
  typedef ImagePairSimilarityMeasureMSD Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Parent class.
  typedef ImagePairSimilarityMeasure Superclass;

  /** Constructor.
   * For reference and model volume, InitDataset is called.
   *@param refVolume The reference (fixed) volume.
   *@param modVolume The model (transformed) volume.
   */
  ImagePairSimilarityMeasureMSD( const UniformVolume::SmartConstPtr& refVolume, const UniformVolume::SmartConstPtr& fltVolume,
				 const Interpolators::InterpolationEnum interpolation = Interpolators::DEFAULT );

  /** Add a pair of values to the metric.
   */
  template<class T> void Increment( const T a, const T b )
  {
    ++this->m_NumberOfSamples;
    this->m_SumOfDifferences -= MathUtil::Square( a - b );
  }

  /** Remove a pair of values from the metric.
   */
  template<class T> void Decrement( const T a, const T b )
  {
    --this->m_NumberOfSamples;
    this->m_SumOfDifferences += MathUtil::Square( a - b );
  }

  /// Reset internal variables for next computation.
  virtual void Reset () 
  {
    this->m_SumOfDifferences = 0;
    this->m_NumberOfSamples = 0;
  }
  
  /// Get the value of the metric.
  Self::ReturnType Get() const 
  {
    return static_cast<Self::ReturnType>( this->m_SumOfDifferences / this->m_NumberOfSamples );
  }

  void Add ( const Self& other )
  {
    this->m_SumOfDifferences += other.m_SumOfDifferences;
    this->m_NumberOfSamples += other.m_NumberOfSamples;
  }

  void Remove ( const Self& other )
  {
    this->m_SumOfDifferences -= other.m_SumOfDifferences;
    this->m_NumberOfSamples -= other.m_NumberOfSamples;
  }
  
private:
  /// this->m_SumOfDifferences of all sample pair differences
  double m_SumOfDifferences;

  /// Number of sample pairs.
  int m_NumberOfSamples;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSimilarityMeasureMSD_h_included_
