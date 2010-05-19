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

#ifdef _MSC_VER
#pragma warning (disable:4521)
#endif
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
    : VoxelMatchingMetric<short,TYPE_SHORT,I>( refVolume, fltVolume ),
    SumJ( NULL ), SumJ2( NULL ), HistogramI( NULL ), 
    SumI( NULL ), SumI2( NULL ), HistogramJ( NULL )
  { 
    NumBinsX = NumBinsY = numBins;

    if ( NumBinsX == CMTK_HISTOGRAM_AUTOBINS )
      NumBinsX = std::max<unsigned>( std::min<unsigned>( refVolume->GetNumberOfPixels(), 128 ), 8 );
    HistogramI = new Histogram<unsigned int>( NumBinsX );

    if ( NumBinsY == CMTK_HISTOGRAM_AUTOBINS )
      NumBinsY = std::max<unsigned>( std::min<unsigned>( fltVolume->GetNumberOfPixels(), 128 ), 8 );
    HistogramJ = new Histogram<unsigned int>( NumBinsY );

    Types::DataItem from, to;
    refVolume->GetData()->GetRange( from, to );
    HistogramI->SetRange( from, to );

    SumJ = Memory::AllocateArray<double>( NumBinsX );
    SumJ2 = Memory::AllocateArray<double>( NumBinsX );

    fltVolume->GetData()->GetStatistics( MuJ, SigmaSqJ );

    fltVolume->GetData()->GetRange( from, to );
    HistogramJ->SetRange( from, to );

    SumI = Memory::AllocateArray<double>( NumBinsY );
    SumI2 = Memory::AllocateArray<double>( NumBinsY );

    refVolume->GetData()->GetStatistics( MuI, SigmaSqI );
  }

  /// Copy constructor.
  VoxelMatchingCorrRatio ( const Self& other )
    : VoxelMatchingMetric<short,TYPE_SHORT,I>( other ) 
  {
    NumBinsX = other.NumBinsX;
    NumBinsY = other.NumBinsY;

    SigmaSqJ = other.SigmaSqJ;
    MuJ = other.MuJ;

    SigmaSqI = other.SigmaSqI;
    MuI = other.MuI;

    HistogramI = new Histogram<unsigned int>( NumBinsX );
    SumJ = Memory::AllocateArray<double>( NumBinsX );
    SumJ2 = Memory::AllocateArray<double>( NumBinsX );

    HistogramJ = new Histogram<unsigned int>( NumBinsY );
    SumI = Memory::AllocateArray<double>( NumBinsY );
    SumI2 = Memory::AllocateArray<double>( NumBinsY );
  }
  
  /// Destructor: free internal data structures.
  ~VoxelMatchingCorrRatio() 
  {
    delete HistogramI;
    delete HistogramJ;
    delete[] SumI;
    delete[] SumI2;
    delete[] SumJ;
    delete[] SumJ2;
  }

  /** Reset computation.
   * Initialize arrays that hold the sums of all floating values and their
   * squares, separated by histogram classes of the reference image.
   */
  void Reset() 
  {
    HistogramI->Reset();
    HistogramJ->Reset();
    memset( SumI, 0, NumBinsY * sizeof( *SumI ) );
    memset( SumI2, 0, NumBinsY * sizeof( *SumI2 ) );
    memset( SumJ, 0, NumBinsX * sizeof( *SumJ ) );
    memset( SumJ2, 0, NumBinsX * sizeof( *SumJ2 ) );
  }

  /** Continue incremental calculation.
   */
  /** Add a pair of values to the metric.
   */
  template<class T> void Increment( const T a, const T b )
  {
    // what's the reference histogram bin?
    size_t bin = HistogramI->ValueToBin( a );
    // count this sample
    HistogramI->Increment( bin );
    // add floating value to sum of values for this class
    SumJ[bin] += b;
    // add squared floating value to sum of squared values for this class
    SumJ2[bin] += b * b;

    // same in reverse
    bin = HistogramJ->ValueToBin( b );
    HistogramJ->Increment( bin );
    SumI[bin] += a;
    SumI2[bin] += a * a;
  }

  /** Remove a pair of values from the metric.
   */
  template<class T> void Decrement( const T a, const T b )
  {
    // what's the reference histogram bin?
    size_t bin = HistogramI->ValueToBin( a );
    // count this sample
    HistogramI->Decrement( bin );
    // add floating value to sum of values for this class
    SumJ[bin] -= b;
    // add squared floating value to sum of squared values for this class
    SumJ2[bin] -= b * b;

    // same in reverse
    bin = HistogramJ->ValueToBin( b );
    HistogramJ->Decrement( bin );
    SumI[bin] -= a;
    SumI2[bin] -= a * a;
  }

  /**
   */
  void AddMetric ( const Self& other )
  {
    HistogramI->AddHistogram( other.HistogramI );
    for ( size_t bin = 0; bin < NumBinsX; ++bin ) 
      {
      SumJ[bin] += other.SumJ[bin];
      SumJ2[bin] += other.SumJ2[bin];
      }
    
    HistogramJ->AddHistogram( other.HistogramJ );
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
    HistogramI->RemoveHistogram( other.HistogramI );
    for ( size_t bin = 0; bin < NumBinsX; ++bin ) {
      SumJ[bin] -= other.SumJ[bin];
      SumJ2[bin] -= other.SumJ2[bin];
    }

    HistogramJ->RemoveHistogram( other.HistogramJ );
    for ( size_t bin = 0; bin < NumBinsY; ++bin ) {
      SumI[bin] -= other.SumI[bin];
      SumI2[bin] -= other.SumI2[bin];
    }
  }

  /// UNDOCUMENTED
  void Copy( const Self& other ) 
  {
    *HistogramI = other.HistogramI;
    for ( size_t bin = 0; bin < NumBinsX; ++bin ) 
      {
      SumJ[bin] = other.SumJ[bin];
      SumJ2[bin] = other.SumJ2[bin];
      }
    SigmaSqJ = other.SigmaSqJ;
    MuJ = other.MuJ;
    
    *HistogramJ = other.HistogramJ;
    for ( size_t bin = 0; bin < NumBinsY; ++bin ) 
      {
      SumI[bin] = other.SumI[bin];
      SumI2[bin] = other.SumI2[bin];
      }
    SigmaSqI = other.SigmaSqI;
    MuI = other.MuI;
  }
  
  /// Return correlation ratio.
  typename Self::ReturnType Get () const;

private:
  /// Number of bins for the X-distribution.
  size_t NumBinsX;
  
  /// Array with sums of all Y-values by X-bins.
  double* SumJ;

  /// Array with sums of squares of all Y-values by X-bins.
  double* SumJ2;

  /// Histogram with counts of all X-values.
  Histogram<unsigned int>* HistogramI;

  // Variance of complete floating image for normalization.
  Types::DataItem SigmaSqJ;

  // Mean of complete floating image for normalization.
  Types::DataItem MuJ;

  /// Number of bins for the Y-distribution.
  size_t NumBinsY;
  
  /// Array with sums of all X-values by Y-bins.
  double* SumI;

  /// Array with sums of squares of all X-values by Y-bins.
  double* SumI2;

  /// Histogram with counts of all X-values.
  Histogram<unsigned int>* HistogramJ;

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
