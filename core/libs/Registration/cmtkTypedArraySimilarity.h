/*
//
//  Copyright 2004-2010 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkTypedArraySimilarity_h_included_
#define __cmtkTypedArraySimilarity_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkTypedArray.h"
#include "Base/cmtkFunctional.h"

#include "Registration/cmtkTypedArraySimilarityMemory.h"

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{
/** Class with operators to compute various pixel similarity measures.
 */
class TypedArraySimilarity 
{
public:
  /// Identifiers for all available metrics.
  typedef enum {
    /// Mutual information.
    MI = 0, 
    /// Normalized mutual information.
    NMI = 1,
    /// Correlation ratio.
    CR = 2,
    /// Cross correlation.
    CC = 3,
    /// Mean squared difference.
    MSD = 4,
    /// Peak signal-to-noise ratio.
    PSNR = 5,
    /// Difference array entropy.
    DAE = 6
  } ID;

  /// Return type for all similarity measures: match cmtk::Functional's
  typedef Functional::ReturnType ReturnType;
  
  /// Compute mutual information between two pixel arrays.
  static ReturnType GetMutualInformation
  ( const TypedArray* array0, const TypedArray* array1,
    TypedArraySimilarityMemory *const memory = NULL );
  
  /// Compute mutual information between two sets of pixel arrays.
  static ReturnType GetMutualInformation( const std::vector<const TypedArray*>& data0, const std::vector<const TypedArray*>& data1, const bool normalized = false );
  
  /// Compute norrmalized mutual information between two sets of pixel arrays.
  static ReturnType GetNormalizedMutualInformation( const std::vector<const TypedArray*>& data0, const std::vector<const TypedArray*>& data1 )
  {
    return GetMutualInformation( data0, data1, true /*normalized*/ );
  }
  
  /** Compute correlation ratio between two pixel arrays.
   * This function is implemented using a 1-D histogram.
   */
  static ReturnType GetCorrelationRatio( const TypedArray* array0, const TypedArray* array1 );

  /// Compute normalized mutual information between two pixel arrays.
  static ReturnType GetNormalizedMutualInformation( const TypedArray* array0, const TypedArray* array1, TypedArraySimilarityMemory *const memory = NULL );

  /// Compute negated (i.e., sign-switched) mean squared pixel difference between two pixel arrays.
  static ReturnType GetMinusMeanSquaredDifference( const TypedArray* array0, const TypedArray* array1 );

  /** Compute Peak-Signal-to-Noise-Ratio.
   *\param data Measured data.
   *\param signal Pure signal without noise.
   */
  static ReturnType GetPeakSignalToNoiseRatio( const TypedArray* data, const TypedArray* signal );

  /// Compute normalized cross correlation between two pixel arrays.
  static ReturnType GetCrossCorrelation( const TypedArray* array0, const TypedArray* array1 );

  /** Compute scaled difference of two images.
   * The values of the second array are scaled with a comon factor so that the
   * entropy of the difference array is minimized. The resulting scale factor
   * is returned via a reference argument.
   */
  static TypedArray::SmartPtr GetDifferenceArray( const TypedArray* array0, const TypedArray* array1, Types::DataItem &scaleFactor );

  /// Compute entropy of difference of two images.
  static ReturnType GetDifferenceArrayEntropy( const TypedArray* array0, const TypedArray* array1, Types::DataItem &scaleFactor );

  /// Check whether two pixel arrays have matching pixel dimensions.
  static bool CheckArrayDimensions( const TypedArray* array0, const TypedArray* array1 );

  /** Compute the optimal scale factor between two images.
   * This implementation uses least squares fitting
   */
  static ReturnType GetOptimalScale( const TypedArray* array0, const TypedArray* array1 );
};

} // namespace

#endif // #ifndef __cmtkTypedArraySimilarity_h_included_

