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

#ifndef __cmtkVoxelMatchingCorrRatio_h_included_
#define __cmtkVoxelMatchingCorrRatio_h_included_

#include <cmtkconfig.h>

#include <cmtkRegistrationJointHistogram.h>
#include <cmtkInterpolator.h>
#include <cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Voxel metric "correlation ratio".
 *\deprecated For future code, use cmtk::ImagePairSimilarityMetricCR instead.
 */
template< Interpolators::InterpolationEnum I = Interpolators::LINEAR >
class VoxelMatchingCorrRatio : 
  public VoxelMatchingMetric<short,TYPE_SHORT,I>
{
public:
  /// This type.
  typedef VoxelMatchingCorrRatio<I> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /** Constructor.
   * The inherited constructor is called to initialize the given datasets.
   * Afterwards, the original (untransformed) probability distribution 
   * functions of model and reference are calculated.
   */
  VoxelMatchingCorrRatio ( const UniformVolume* refVolume, const UniformVolume* fltVolume, const unsigned int numBins = CMTK_HISTOGRAM_AUTOBINS )
    : VoxelMatchingMetric<short,TYPE_SHORT,I>( refVolume, fltVolume )
  { 
    NumBinsX = NumBinsY = numBins;

    if ( NumBinsX == CMTK_HISTOGRAM_AUTOBINS )
      NumBinsX = std::max<unsigned>( std::min<unsigned>( refVolume->GetNumberOfPixels(), 128 ), 8 );
    HistogramI.Resize( NumBinsX );

    if ( NumBinsY == CMTK_HISTOGRAM_AUTOBINS )
      NumBinsY = std::max<unsigned>( std::min<unsigned>( fltVolume->GetNumberOfPixels(), 128 ), 8 );
    HistogramJ.Resize( NumBinsY );

    HistogramI.SetRange( refVolume->GetData()->GetRange() );

    SumJ.resize( NumBinsX );
    SumJ2.resize( NumBinsX );
    
    fltVolume->GetData()->GetStatistics( MuJ, SigmaSqJ );

    HistogramJ.SetRange( fltVolume->GetData()->GetRange() );

    SumI.resize( NumBinsY );
    SumI2.resize( NumBinsY );
    
    refVolume->GetData()->GetStatistics( MuI, SigmaSqI );
  }

  /** Reset computation.
   * Initialize arrays that hold the sums of all floating values and their
   * squares, separated by histogram classes of the reference image.
   */
  void Reset() 
  {
    HistogramI.Reset();
    HistogramJ.Reset();
    std::fill( SumI.begin(), SumI.end(), 0 );
    std::fill( SumJ.begin(), SumJ.end(), 0 );
    std::fill( SumI2.begin(), SumI2.end(), 0 );
    std::fill( SumJ2.begin(), SumJ2.end(), 0 );
  }

  /** Continue incremental calculation.
   */
  /** Add a pair of values to the metric.
   */
  template<class T> void Increment( const T a, const T b )
  {
    // what's the reference histogram bin?
    size_t bin = HistogramI.ValueToBin( a );
    // count this sample
    HistogramI.Increment( bin );
    // add floating value to sum of values for this class
    SumJ[bin] += b;
    // add squared floating value to sum of squared values for this class
    SumJ2[bin] += b * b;

    // same in reverse
    bin = HistogramJ.ValueToBin( b );
    HistogramJ.Increment( bin );
    SumI[bin] += a;
    SumI2[bin] += a * a;
  }

  /** Remove a pair of values from the metric.
   */
  template<class T> void Decrement( const T a, const T b )
  {
    // what's the reference histogram bin?
    size_t bin = HistogramI.ValueToBin( a );
    // count this sample
    HistogramI.Decrement( bin );
    // add floating value to sum of values for this class
    SumJ[bin] -= b;
    // add squared floating value to sum of squared values for this class
    SumJ2[bin] -= b * b;

    // same in reverse
    bin = HistogramJ.ValueToBin( b );
    HistogramJ.Decrement( bin );
    SumI[bin] -= a;
    SumI2[bin] -= a * a;
  }

  /**
   */
  void AddMetric ( const Self& other )
  {
    HistogramI.AddHistogram( other.HistogramI );
    for ( size_t bin = 0; bin < NumBinsX; ++bin ) 
      {
      SumJ[bin] += other.SumJ[bin];
      SumJ2[bin] += other.SumJ2[bin];
      }
    
    HistogramJ.AddHistogram( other.HistogramJ );
    for ( size_t bin = 0; bin < NumBinsY; ++bin ) 
      {
      SumI[bin] += other.SumI[bin];
      SumI2[bin] += other.SumI2[bin];
      }
  }
  
  /**
   */
  void RemoveMetric ( const Self& other )
  {
    HistogramI.RemoveHistogram( other.HistogramI );
    for ( size_t bin = 0; bin < NumBinsX; ++bin ) {
      SumJ[bin] -= other.SumJ[bin];
      SumJ2[bin] -= other.SumJ2[bin];
    }

    HistogramJ.RemoveHistogram( other.HistogramJ );
    for ( size_t bin = 0; bin < NumBinsY; ++bin ) {
      SumI[bin] -= other.SumI[bin];
      SumI2[bin] -= other.SumI2[bin];
    }
  }

  /// Return correlation ratio.
  typename Self::ReturnType Get () const;

private:
  /// Number of bins for the X-distribution.
  size_t NumBinsX;
  
  /// Array with sums of all Y-values by X-bins.
  std::vector<double> SumJ;

  /// Array with sums of squares of all Y-values by X-bins.
  std::vector<double> SumJ2;

  /// Histogram with counts of all X-values.
  Histogram<unsigned int> HistogramI;

  // Variance of complete floating image for normalization.
  Types::DataItem SigmaSqJ;

  // Mean of complete floating image for normalization.
  Types::DataItem MuJ;

  /// Number of bins for the Y-distribution.
  size_t NumBinsY;
  
  /// Array with sums of all X-values by Y-bins.
  std::vector<double> SumI;

  /// Array with sums of squares of all X-values by Y-bins.
  std::vector<double> SumI2;

  /// Histogram with counts of all X-values.
  Histogram<unsigned int> HistogramJ;

  // Variance of complete floating image for normalization.
  Types::DataItem SigmaSqI;

  // Mean of complete floating image for normalization.
  Types::DataItem MuI;
};

/// Correlation ratio with trilinear interpolation.
typedef VoxelMatchingCorrRatio<Interpolators::LINEAR> VoxelMatchingCorrRatio_Trilinear;

/// Correlation ratio with nearest-neighbor interpolation.
typedef VoxelMatchingCorrRatio<Interpolators::NEAREST_NEIGHBOR> VoxelMatchingCorrRatio_NearestNeighbor;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVoxelMatchingCorrRatio_h_included_
