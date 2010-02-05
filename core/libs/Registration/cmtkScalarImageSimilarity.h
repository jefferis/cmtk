/*
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkScalarImageSimilarity_h_included_
#define __cmtkScalarImageSimilarity_h_included_

#include <cmtkconfig.h>

#include <cmtkScalarImage.h>
#include <cmtkTypedArraySimilarity.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Memory for image similarity computation.
class ScalarImageSimilarityMemory :
  /// This is basically the same as the typed array similarity memory.
  public TypedArraySimilarityMemory
{
};

/** Class with operators to compute various 2-D image similarity measures.
 * Most functions in this class implement operators investigated by Penney
 * et al. [Penney G, et al., Similarity Measures for Use in 2D-3D Medical Image
 * Registration, IEEE Trans Med Imaging 17(4):586-595, 1998 (August)] for use
 * in registration of 2-D fluoroscopic images to 3-D CT volumes.
 *@see TypedArraySimilarity for implementation of all operators that do not
 * require 2-D arrangement of pixels but work on the sequential pixel arrays.
 */
class ScalarImageSimilarity :
  /// Inherit similarity functions for unstructured arrays.
  public TypedArraySimilarity
{
public:
  /// Identifiers for all available metrics.
  typedef enum {
    /// Mutual information.
    MI = 0, 
    /// Normalized mutual information.
    NMI = 1,
    /// Regional mutual information.
    RMI = 2,
    /// Regional normalized mutual information.
    RNMI = 3,
    /// Correlation ratio.
    CR = 4,
    /// Cross correlation.
    CC = 5,
    /// Mean squared difference.
    MSD = 6,
    /// Peak signal-to-noise ratio.
    PSNR = 7,
    /// Difference array entropy.
    DAE = 8,
    /// Gradient correlation.
    GradientCorrelation = 9,
    /// Pattern intensity.
    PatternIntensity = 10
  } ID;
  
  /** Compute mutual information between two images.
   * This is the standard mutual information similarity measure. See Penney et 
   * al. III-C.
   */
  static ReturnType GetMutualInformation( const ScalarImage* image0, const ScalarImage* image1, ScalarImageSimilarityMemory *const memory = NULL );

  /** Compute normalized mutual information between two images.
   * This is the normalized mutual information measure as described by
   * Studholme et al. [Studholme, Hill, Hawkes, Pattern Analysis, 1999].
   */
  static ReturnType GetNormalizedMutualInformation( const ScalarImage* image0, const ScalarImage* image1, ScalarImageSimilarityMemory *const memory = NULL );

  /// Compute mutual information between two sets of pixel arrays.
  static ReturnType GetMutualInformation( std::vector<const ScalarImage*>& array0, std::vector<const ScalarImage*>& array1, const bool normalized = false )
  {
    std::vector<const TypedArray*> data0(array0.size());
    std::vector<const TypedArray*> data1(array1.size());
    for ( size_t i = 0; i < array0.size(); ++i )
      {
      data0[i] = array0[i]->GetPixelData();
      }
    for ( size_t i = 0; i < array1.size(); ++i )
      {
      data1[i] = array1[i]->GetPixelData();
      } 
    return TypedArraySimilarity::GetMutualInformation( data0, data1, normalized );
  }
  
  /// Compute norrmalized mutual information between two sets of pixel arrays.
  static ReturnType GetNormalizedMutualInformation( std::vector<const ScalarImage*>& array0, std::vector<const ScalarImage*>& array1 )
  {
    return GetMutualInformation( array0, array1, true /*normalized*/ );
  }
  
  /** Compute regional mutual information between two images.
   */
  static ReturnType GetRegionalMutualInformation( const ScalarImage* image0, const ScalarImage* image1, const int radius = 2 );

  /** Compute negated (i.e., sign-switched) mean squared pixel difference between two images.
   */
  static ReturnType GetMinusMeanSquaredDifference( const ScalarImage* image0, const ScalarImage* image1 );

  /** Compute normalized cross correlation between two images.
   * See Penney et al. III-A.
   */
  static ReturnType GetCrossCorrelation( const ScalarImage* image0, const ScalarImage* image1 );

  /** Compute gradient correlation between two images.
   * See Penney et al. III-D.
   */
  static ReturnType GetGradientCorrelation( const ScalarImage* image0, const ScalarImage* image1 );

  /** Compute gradient difference between two images.
   * See Penney et al. III-F.
   */
  static ReturnType GetGradientDifference( const ScalarImage* image0, const ScalarImage* image1, const ReturnType Ax = 1, const ReturnType Ay = 1 );

  /** Compute pattern intensity between two images.
   * See Penney et al. III-E.
   */
  static ReturnType GetPatternIntensity( const ScalarImage* image0, const ScalarImage* image1, const ReturnType sigma = 10, const int radius = 3 );

  /** Compute entropy of difference of two images.
   * See Penney et al. III-B.
   */
  static ReturnType GetDifferenceImageEntropy( const ScalarImage* image0, const ScalarImage* image1 );

  /** Compute entropy of difference of two images.
   * See Penney et al. III-B. This function returns the actual scaling factor
   * used for intensity normalization via a reference parameter.
   */
  static ReturnType GetDifferenceImageEntropy( const ScalarImage* image0, const ScalarImage* image1, Types::DataItem &scaleFactor );


  /** Compute correlation ratio between two pixel arrays.
   * This function is implemented using a 1-D histogram.
   */
  static ReturnType GetCorrelationRatio( const ScalarImage* image0, const ScalarImage* image1 );

  /// Check whether two images are valid and have matching pixel dimensions.
  static bool CheckImageDimensions( const ScalarImage* image0, const ScalarImage* image1 );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkScalarImageSimilarity_h_included_
