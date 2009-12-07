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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkImagePairSimilarityMeasureNCC_h_included_
#define __cmtkImagePairSimilarityMeasureNCC_h_included_

#include <cmtkconfig.h>

#include <cmtkImagePairSimilarityMeasure.h>
#include <cmtkUniformVolume.h>
#include <cmtkTypedArray.h>

#include <cmtkMathUtil.h>

#include <cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

#ifdef _MSC_VER
#pragma warning (disable:4521)
#endif
/** Normalized Cross Correlation Metric.
 */
class ImagePairSimilarityMeasureNCC :
  /// Inherit generic pairwise similarity measure class
  public ImagePairSimilarityMeasure
{
public:
  /// This type.
  typedef ImagePairSimilarityMeasureNCC Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /** Default constructor.
   */
  ImagePairSimilarityMeasureNCC() {};

  /** Constructor.
   * For reference and floating volume, InitDataset is called.
   *@param refVolume The reference (fixed) volume.
   *@param fltVolume The floating (transformed) volume.
   */
  ImagePairSimilarityMeasureNCC ( const UniformVolume::SmartPtr& refVolume, const UniformVolume::SmartPtr& fltVolume,
				  const Interpolators::InterpolationEnum interpolation = Interpolators::DEFAULT );

  /// Copy constructor.
  ImagePairSimilarityMeasureNCC ( Self& other, const bool copyData = false );
  
  /// Constant copy constructor.
  ImagePairSimilarityMeasureNCC ( const Self& other );
  
  /// Copy.
  void Copy ( const ImagePairSimilarityMeasureNCC& other, const bool copyData = false ) 
  {
    this->CopyUnsafe( other, copyData );
  }
  
  /// Copy values, but do not create or check data structures.
  void CopyUnsafe ( const ImagePairSimilarityMeasureNCC& other, const bool copyData = false );

  /** Add a pair of values to the metric.
   */
  template<class T> void Increment( const T a, const T b )
  {
    ++Samples;
    SumX += a;
    SumY += b;
    SumSqX += a * a;
    SumSqY += b * b;
    SumXY += a * b;
  }

  /** Remove a pair of values from the metric.
   */
  template<class T> void Decrement( const T a, const T b )
  {
    --Samples;
    SumX -= a;
    SumY -= b;
    SumSqX -= a * a;
    SumSqY -= b * b;
    SumXY -= a * b;
  }
  
  /// Start with a new computation.
  void Reset () 
  {
    SumX = SumY = SumSqX = SumSqY = SumXY = 0;
    Samples = 0;
  }
  
  /// Compute cross correlation.
  Self::ReturnType Get() const;

  void Add ( const Self& other )
  {
    SumX += other.SumX;
    SumY += other.SumY;
    SumXY += other.SumXY;
    SumSqX += other.SumSqX;
    SumSqY += other.SumSqY;
    Samples += other.Samples;
  }

  void Remove ( const Self& other )
  {
    SumX -= other.SumX;
    SumY -= other.SumY;
    SumXY -= other.SumXY;
    SumSqX -= other.SumSqX;
    SumSqY -= other.SumSqY;
    Samples -= other.Samples;
    assert( Samples >= 0 );
  }

private:
  /// Sum over all samples in X distribution.
  double SumX;

  /// Sum over all samples in Y distribution.
  double SumY;

  /// Sum over products of corresponding samples in X and Y distribution.
  double SumXY;

  /// Sum over all squared samples in X distribution.
  double SumSqX;

  /// Sum over all squared samples in Y distribution.
  double SumSqY;

  /// Number of samples.
  size_t Samples;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImagePairSimilarityMeasureNCC_h_included_
