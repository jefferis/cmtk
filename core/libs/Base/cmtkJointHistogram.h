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

#ifndef __cmtkJointHistogram_h_included_
#define __cmtkJointHistogram_h_included_

#include <cmtkconfig.h>

#include <cmtkJointHistogramBase.h>

#include <cmtkUniformVolume.h>
#include <cmtkTypedArray.h>
#include <cmtkMathUtil.h>

#include <cmtkHistogram.h>

#include <stdio.h>

#include <cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Joint histogram of two distributions.
 *@param T Template parameter: the type of the histogram bins. Can be integral,
 * or float in case of fractional bins.
 */
template<class T>
class JointHistogram : 
  /// Inherit from non-templated base class.
  public JointHistogramBase 
{
protected:
  /// Number of reference data bins.
  size_t NumBinsX;

  /// Width of reference data bins.
  Types::DataItem BinWidthX;

  /// Lower value bound of reference data bins.
  Types::DataItem BinOffsetX;

  /// Number of model data bins.
  size_t NumBinsY;

  /// Width of model data bins.
  Types::DataItem BinWidthY;

  /// Lower value bound of model data bins.
  Types::DataItem BinOffsetY;

  /// Array of cross-modality (joint) bins.
  T* JointBins;
  
public:
  /// This class.
  typedef JointHistogram<T> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /** Default constructor.
   */
  JointHistogram () 
  {
    RealNumBinsX = RealNumBinsY = NumBinsX = NumBinsY = 0; 
    BinWidthX = BinWidthY = 1.0;
    BinOffsetX = BinOffsetY = 0.0;
    this->JointBins = NULL; 
  }

  /** Constructor.
   */
  JointHistogram( const size_t numBinsX, const size_t numBinsY, const bool reset = true ) 
  {
    this->JointBins = NULL;
    RealNumBinsX = RealNumBinsY = NumBinsX = NumBinsY = 0; 
    BinWidthX = BinWidthY = 1.0;
    BinOffsetX = BinOffsetY = 0.0;
    this->SetNumBins( numBinsX, numBinsY, reset );
  }

  /// Copy constructor.
  JointHistogram( const Self& other, const bool copyData = true ) 
  {
    this->JointBins = NULL;
    RealNumBinsX = RealNumBinsY = NumBinsX = NumBinsY = 0; 
    BinWidthX = other.BinWidthX;
    BinWidthY = other.BinWidthY;
    BinOffsetX = other.BinOffsetX;
    BinOffsetY = other.BinOffsetY;

    this->SetNumBins( other.NumBinsX, other.NumBinsY, !copyData /*reset*/ );
    if ( other.JointBins && copyData )
      memcpy( this->JointBins, other.JointBins, RealNumBinsX * RealNumBinsY * sizeof( T ) );
  }

  /** Destructor.
   * All bin arrays and the precomputed data bin index arrays are
   * de-allocated.
   */
  ~JointHistogram () 
  {  
    if ( this->JointBins ) 
      delete[] this->JointBins; 
  }

  /// Resize and allocate histogram bins.
  void Resize( const size_t numberOfBinsX, const size_t numberOfBinsY, const bool reset = true )
  {
    if ( (numberOfBinsX != this->NumBinsX) || (numberOfBinsY != this->NumBinsY) )
      {
      this->SetNumBins( numberOfBinsX, numberOfBinsY );
      }
    if ( reset ) 
      this->Reset();
  }

  /** Make an identical copy of this object.
   *@param copyData If this is non-zero (default), then the current values in
   * the histogram bins are also duplicated. If zero, only the histogram
   * structure is duplicated.
   *@return The newly created copy of this object.
   */
  Self* Clone ( const bool copyData = true ) const 
  {
    Self *clone = new Self( NumBinsX, NumBinsY, false );
    
    if ( copyData )
      memcpy( clone->JointBins, this->JointBins, RealNumBinsX * RealNumBinsY * sizeof( T ) );
    else
      clone->Reset();
    
    return clone;
  }

  /// Copy another histogram after matching target size.
  void Copy ( const Self& other ) 
  {
    this->Resize( other.NumBinsX, other.NumBinsY, false /*reset*/ );
    this->CopyUnsafe( other );
  }
  
  /// Copy another histogram without range checking.
  void CopyUnsafe ( const Self& other ) 
  {
    const size_t realNumBinsXY = this->RealNumBinsX * this->RealNumBinsY;
    for ( size_t idx = 0; idx < realNumBinsXY; ++idx )
      this->JointBins[idx] = other.JointBins[idx];
//    memcpy( this->JointBins, other.JointBins, RealNumBinsX * RealNumBinsY * sizeof( T ) );
  }

  /// Set histogram values from existing data array.
  void SetBins( const T* bins ) 
  {
    size_t fromOffset = 0, toOffset = 0;
    for ( size_t y=0; y < NumBinsY; ++y ) 
      {
      for ( size_t x=0; x < NumBinsX; ++x, ++fromOffset, ++toOffset )
	this->JointBins[ toOffset ] = bins[ fromOffset ];
      this->JointBins[ toOffset++ ] = 0;
      }
    for ( size_t x=0; x < RealNumBinsX; ++x, ++toOffset ) 
      this->JointBins[ toOffset ] = 0;
  }

  /// Get histogram values and put them into data array.
  T* GetBins() const 
  {
    T *bins = Memory::AllocateArray<T>( NumBinsX * NumBinsY );
    
    size_t fromOffset = 0, toOffset = 0;
    for ( size_t y=0; y < NumBinsY; ++y ) 
      {
      for ( size_t x=0; x < NumBinsX; ++x, ++fromOffset, ++toOffset )
	bins[ toOffset ] = this->JointBins[ fromOffset ];
      ++fromOffset;
      }
    return bins;
  }
  
  /// Get constant pointer to internal histogram values array.
  const T* GetRawBins() const 
  {
    return this->JointBins;
  }

  /** Set number of bins in x-direction.
   * This will re-allocate the histogram storage and clear all entries.
   */
  void SetNumBinsX( const size_t numBinsX ) 
  {
    this->SetNumBins( numBinsX, NumBinsY );
  }

  /** Set number of bins in y-direction.
   * This will re-allocate the histogram storage and clear all entries.
   */
  void SetNumBinsY( const size_t numBinsY ) 
  {
    this->SetNumBins( NumBinsX, numBinsY );
  }
  
  /** Set number of bins in x- and y-direction.
   * This will re-allocate the histogram storage and clear all entries.
   */
  void SetNumBins
  ( const size_t numBinsX, const size_t numBinsY, const bool reset = true ) 
  {
    if ( this->JointBins ) delete[] this->JointBins;
    
    this->NumBinsX = numBinsX;
    this->NumBinsY = numBinsY;
    
    this->RealNumBinsX = NumBinsX + 1;
    this->RealNumBinsY = NumBinsY + 1;
    
    this->JointBins = Memory::AllocateArray<T>( this->RealNumBinsX * this->RealNumBinsY );
    if ( reset ) this->Reset();
  }

  /// Return number of histogram bins in X direction.
  size_t GetNumBinsX() const 
  { 
    return this->NumBinsX;
  }

  /// Return number of histogram bins in Y direction.
  size_t GetNumBinsY() const
  { 
    return this->NumBinsY;
  }

  /** Set value range of the X distribution.
   */
  void SetRangeX( const Types::DataItem from, const Types::DataItem to ) 
  {
    this->BinOffsetX = from;
    this->BinWidthX = ( to - from ) / ( this->NumBinsX - 1 );
  }

  /** Set value range of the Y distribution.
   */
  void SetRangeY( const Types::DataItem from, const Types::DataItem to ) 
  {
    this->BinOffsetY = from;
    this->BinWidthY = ( to - from ) / ( this->NumBinsY - 1 );
  }
  
  /** Set value range of the X distribution where upper and lower bound are centered in first and last histogram bins.
   */
  void SetRangeCenteredX( const Types::DataItem from, const Types::DataItem to ) 
  {
    this->BinWidthX = ( to - from ) / (this->NumBinsX - 1);
    this->BinOffsetX = - this->BinWidthX / 2;
  }

  /** Set value range of the Y distribution where upper and lower bound are centered in first and last histogram bins.
   */
  void SetRangeCenteredY( const Types::DataItem from, const Types::DataItem to ) 
  {
    this->BinWidthY = ( to - from ) / (this->NumBinsY - 1);
    this->BinOffsetY = - this->BinWidthY / 2;
  }
  
  /** Get value range of the X distribution.
   */
  void GetRangeX( Types::DataItem& from, Types::DataItem& to ) const 
  {
    from = this->BinOffsetX;
    to = from + this->BinWidthX * ( this->NumBinsX - 1);
  }
  
  /** Get value range of the Y distribution.
   */
  void GetRangeY( Types::DataItem& from, Types::DataItem& to ) const 
  {
    from = this->BinOffsetY;
    to = from + this->BinWidthY * ( this->NumBinsY - 1 );
  }
  
  /** Reset computation.
   * This function has to be called before any other computation made with an
   * object of this class. All bin counters are reset to zero, therefore
   * Reset() must also be called before any new computation performed using an
   * already previously used object.
   */
  void Reset () 
  {
    memset( this->JointBins, 0, RealNumBinsX * RealNumBinsY * sizeof( T ) );
  }
  
  /** Return bin corresponding to a certain value of the X distribution.
   *@param value A value from the X distribution.
   *@return The index of the X bin corresponding to the given value.
   */
  size_t ValueToBinX ( const Types::DataItem value ) const 
  {
    size_t binIndex = 
      static_cast<size_t>( (value - BinOffsetX) / BinWidthX );
    return std::max<size_t>( 0, std::min<size_t>( NumBinsX-1, binIndex ) );
  }
  
  /** Return bin corresponding to a certain value of the Y distribution.
   *@param value A value from the Y distribution.
   *@return The index of the Y bin corresponding to the given value.
   */
  size_t ValueToBinY ( const Types::DataItem value ) const 
  {
    size_t binIndex = 
      static_cast<int>( (value - BinOffsetY) / BinWidthY );
    return std::max<size_t>( 0, std::min<size_t>( NumBinsY-1, binIndex ) );
  }

  /** Return center of values represented by a certain X distribution bin.
   *@param bin Index of a bin from the X distribution.
   *@return Average of upper and lower margin values of the given bin.
   */
  Types::DataItem BinToValueX ( const size_t bin ) const 
  {
    return BinOffsetX + (bin+0.5) * BinWidthX;
  }

  /** Return center of values represented by a certain Y distribution bin.
   *@param bin Index of a bin from the Y distribution.
   *@return Average of upper and lower margin values of the given bin.
   */
  Types::DataItem BinToValueY ( const size_t bin ) const 
  {
    return BinOffsetY + (bin+0.5) * BinWidthY;
  }

  /** Return projection of 2D distribution to X.
   * This function can be used to reconstruct the marginal distribution X from
   * the 2D histogram without explicitly stored 1D distributions.
   *@param indexX A bin index for the X distribution.
   *@return Projection of 2D histogram onto X distribution. This is the sum
   * of all bins in Y direction for the given index in X direction.
   */
  T ProjectToX ( const size_t indexX ) const 
  {
    T project = 0;
    for ( size_t j=0; j < NumBinsY; ++j )
      project += this->JointBins[indexX + j * RealNumBinsX];
    
    return project;
  }
  
  /** Return projection of 2D distribution to Y.
   * This function can be used to reconstruct the marginal distribution Y from
   * the 2D histogram without explicitly stored 1D distributions.
   *@param indexY A bin index for the Y distribution.
   *@return Projection of 2D histogram onto X distribution. This is the sum
   * of all bins in X direction for the given index in X direction.
   */
  T ProjectToY ( const size_t indexY ) const 
  {
    T project = 0;
    size_t offset = indexY * RealNumBinsX;
    for ( size_t i = 0; i < NumBinsX; ++i )
      project += this->JointBins[i + offset];
    
    return project;
  }

  /// Get marginal X distribution as 1D histogram.
  Histogram<T>* GetMarginalX() const;

  /// Get marginal Y distribution as 1D histogram.
  Histogram<T>* GetMarginalY() const;

  /** Return a bin value.
   */
  T GetBin( const size_t indexX, const size_t indexY ) const 
  {
    return this->JointBins[ indexX + RealNumBinsX * indexY ];
  }
  
  /** Return total number of samples stored in the histogram.
   */
  T SampleCount () const 
  {
    T sampleCount = 0;
    
    size_t idx = 0;
    for ( size_t i=0; i<NumBinsY; ++i ) 
      {
      for ( size_t j=0; j<NumBinsX; ++j, ++idx )
	sampleCount += this->JointBins[idx];
      
      ++idx; // Skip extra bin at the end of each row.
      }
    
    return sampleCount;
  }

  /** Return maximum number of samples stored in any bin.
   */
  T GetMaximumBinValue () const 
  {
    T maximum = 0;
    
    size_t idx = 0;
    for ( size_t i=0; i<NumBinsY; ++i ) 
      {
      for ( size_t j=0; j<NumBinsX; ++j, ++idx )
	maximum = std::max( maximum, this->JointBins[idx] );
      
      /// Skip extra bin at the end of each row.
	++idx;
      }
    
    return maximum;
  }
  
  /** Compute marginal entropies.
   * From the bin counts, the marginal entropies of both, reference and
   * model data are estimated.
   *@param HX Upon return, this reference holds the estimated marginal entropy
   * of the X random variable, i.e. the reference image.
   *@param HY Upon return, this reference holds the estimated marginal entropy
   * of the Y random variable, i.e. the model image.
   */
  void GetMarginalEntropies ( double& HX, double& HY ) const;

  /** Compute joint entropy.
   * From the joint bins, an estimate to the joint entropy of both, 
   * reference and model image is computed.
   *@param HXY Upon return, this reference holds the estimated joint entropy.
   */
  void GetJointEntropy ( double& HXY ) const
  {
    HXY = this->GetJointEntropy();
  }
  
  /** Compute joint entropy.
   * From the joint bins, an estimate to the joint entropy of both, 
   * reference and model image is computed.
   *@param HXY Upon return, this reference holds the estimated joint entropy.
   */
  double GetJointEntropy() const;
  
  /** Increment the value of a histogram bin by 1.
   * The histogram field to increment is identified directly by its index;
   * no value-rescaling is done internally.
   *@param sampleX Index of histogram field in x-direction.
   *@param sampleY Index of histogram field in y-direction.
   */
  void Increment ( const size_t sampleX, const size_t sampleY ) 
  {
    ++this->JointBins[sampleX+sampleY*RealNumBinsX];
  }

  /** Increment the value of a histogram bin by an arbitrary value.
   * The histogram field to increment is identified directly by its index;
   * no value-rescaling is done internally.
   *@param sampleX Index of histogram field in x-direction.
   *@param sampleY Index of histogram field in y-direction.
   *@param weight Value to increment the given histogram bin by.
   */
  void Increment 
  ( const size_t sampleX, const size_t sampleY, const T weight ) 
  {
    this->JointBins[ sampleX + sampleY * RealNumBinsX ] += weight;
  }

  /** Decrement the value of a histogram bin by 1.
   * The histogram field to decrement is identified directly by its index;
   * no value-rescaling is done internally. Make sure that a value has actually
   * been added to this bin before - otherwise, the next entropy computation my
   * give some unexpected results.
   *@param sampleX Index of histogram field in x-direction.
   *@param sampleY Index of histogram field in y-direction.
   */
  void Decrement ( const size_t sampleX, const size_t sampleY ) 
  {
    --this->JointBins[sampleX+sampleY*RealNumBinsX];
  }
  
  void Decrement
  ( const size_t sampleX, const size_t sampleY, const Types::DataItem weight ) 
  {
    this->JointBins[ sampleX + sampleY * RealNumBinsX ] -= static_cast<T>( weight );
  }
  
  /** Add values from another histogram.
   * Adding is done by corresponding bins. The caller has to make sure that
   * both histograms actually have the same number and arrangement of bins.
   * It is also a good idea to ensure that the data ranges of these bins are
   * the same in both objects. Both can be guaranteed if one histogram was
   * created from the other by a call to Clone() for example.
   *@param other A pointer to the other histogram. Its bin values are added to
   * this object's bins.
   */
  void AddJointHistogram ( const Self& other ) 
  {
    const size_t realNumBinsXY = this->RealNumBinsX * this->RealNumBinsY;
    for ( size_t idx = 0; idx < realNumBinsXY; ++idx )
      this->JointBins[idx] += other.JointBins[idx];
  }
  
  /** Add values from another 1D histogram in X-direction.
   * Adding is done by corresponding bins. The caller has to make sure that
   * both histograms actually have the same number and arrangement of bins.
   * It is also a good idea to ensure that the data ranges of these bins are
   * the same in both objects.
   *@param other A pointer to the other histogram. Its bin values are added to
   * this object's bins.
   *@param sampleY Index of the histogram row to which the 1D histogram is to
   * be added.
   */
  void AddHistogramRow( const Histogram<T>& other, const size_t sampleY, const float weight = 1 ) 
  {
    size_t idx = RealNumBinsX * sampleY;
    for ( size_t i = 0; i<NumBinsX; ++i, ++idx )
      {
      this->JointBins[idx] += static_cast<T>( weight * other.GetBin( i ) );
      }
  }
  
  /** Add values from another 1D histogram in Y-direction.
   * Adding is done by corresponding bins. The caller has to make sure that
   * both histograms actually have the same number and arrangement of bins.
   * It is also a good idea to ensure that the data ranges of these bins are
   * the same in both objects.
   *@param other A pointer to the other histogram. Its bin values are added to
   * this object's bins.
   *@param sampleX Index of the histogram column to which the 1D histogram is
   * to be added.
   */
  void AddHistogramColumn( const size_t sampleX, const Histogram<T>& other, const float weight = 1 ) 
  {
    size_t idx = sampleX;
    for ( size_t j = 0; j<NumBinsY; ++j, idx += RealNumBinsX )
      this->JointBins[idx] += static_cast<T>( weight * other.GetBin( j ) );
  }
  
  /** Subtract bin values from another histogram.
   * Subtraction is done by corresponding bins. The caller has to make sure 
   * that both histograms actually have the same number and arrangement of 
   * bins. It is also a good idea to ensure that the data ranges of these bins
   * are the same in both objects. Both can be guaranteed if one histogram was
   * created from the other by a call to Clone() for example.
   *@param other A pointer to the other histogram. Its bin values are
   * subtracted this object's bins.
   */
  void RemoveJointHistogram ( const Self& other ) 
  {
    const size_t realNumBinsXY = this->RealNumBinsX * this->RealNumBinsY;
    for ( size_t idx = 0; idx < realNumBinsXY; ++idx )
      {
      this->JointBins[idx] -= other.JointBins[idx];
      }
  }
  
  /** Normalize histogram values over X dimension.
   *@param normalizeTo All histogram bins in every row of the histogram are
   * scaled by a common factor so that their sum matches the value of this 
   * parameter.
   */
  void NormalizeOverX( const double normalizeTo = 1.0 ) 
  {
    for ( size_t j = 0; j < NumBinsY; ++j ) 
      {
      const T project = this->ProjectToY( j );
      if ( project > 0 )
	{
	const double factor = normalizeTo / project;
	for ( size_t i = 0; i < NumBinsX; ++i )
	  this->JointBins[ i + RealNumBinsX * j ] = static_cast<T>( this->JointBins[ i + RealNumBinsX * j ] * factor );
	}
      }
  }
  
  /** Normalize histogram values over Y dimension.
   *@param normalizeTo All histogram bins in every column of the histogram are
   * scaled by a common factor so that their sum matches the value of this 
   * parameter.
   */
  void NormalizeOverY( const double normalizeTo = 1.0 ) 
  {
    for ( size_t i = 0; i < NumBinsX; ++i ) 
      {
      const T project = this->ProjectToX( i );
      if ( project > 0 ) 
	{
	const double factor = normalizeTo / project;
	for ( size_t j = 0; j < NumBinsY; ++j )
	  this->JointBins[ i + RealNumBinsX * j ] =
	    static_cast<T>( this->JointBins[ i + RealNumBinsX * j ] * factor );
	}
      }
  }
  
  /* Return the index of the bin with the maximum value for one row.
   *@param j Index of the row.
   *@return The index of the bin with the maximum value in row j.
   */
  size_t GetMaximumBinIndexOverX( const size_t j ) const 
  {
    size_t offset = j * RealNumBinsX;
    
    size_t maxIndex = 0;
    T maxValue = this->JointBins[ offset ];
    
    for ( size_t i = 1; i < NumBinsX; ++i ) 
      {
      offset++;
      if ( this->JointBins[ offset ] > maxValue ) 
	{
	maxValue = this->JointBins[ offset ];
	maxIndex = i;
	}
      }
    return maxIndex;
  }
  
  /* Return the index of the bin with the maximum value for one column.
   *@param j Index of the column.
   *@return The index of the bin with the maximum value in column j.
   */
  size_t GetMaximumBinIndexOverY( const size_t i ) const 
  {
    size_t offset = i;
    
    size_t maxIndex = 0;
    T maxValue = this->JointBins[ offset ];
    
    for ( size_t j = 1; j < NumBinsY; ++j ) 
      {
      offset += RealNumBinsX;
      if ( this->JointBins[ offset ] > maxValue ) 
	{
	maxValue = this->JointBins[ offset ];
	maxIndex = j;
	}
      }
    return maxIndex;
  }

  /** Compute the (normalized) mutual information of this histogram.
   *@param normalized If this parameter is true, then the "normalized" version
   * of Mutual Information as introduced by Studholme et al. [Pattern Analysis
   * 1999] is computed. If this flag is false, the original formula of Viola
   * [IEEE Trans Med Imaging 1997] and Maes is used instead.
   */
  double GetMutualInformation( const bool normalized = false ) const 
  {
    double hX, hY, hXY;
    this->GetMarginalEntropies( hX, hY );
    this->GetJointEntropy( hXY );
    if ( hXY > 0 )
      if ( normalized ) 
	return (hX + hY) / hXY;
      else
	return (hX + hY) - hXY;
    else
      return 0;
  }

  /// Compute the Correlation Ratio corresponding to this histogram.
  double GetCorrelationRatio() const 
  {
    T sampleCount = this->SampleCount();
    if ( ! sampleCount ) return 0;
    
    double sigSquare = 0, m = 0;
    for ( size_t j = 0; j < NumBinsY; ++j ) 
      {
      sigSquare += ( MathUtil::Square(j) * ProjectToY(j) );
      m += (j * ProjectToY(j));
      }
    
    double invSampleCount = 1.0 / sampleCount;
    
    m *= invSampleCount;
    (sigSquare *= invSampleCount) -= MathUtil::Square(m);
    
    double eta = 0;
    for ( size_t i = 0; i < NumBinsX; ++i ) 
      {
      if ( ProjectToX(i) > 0 ) 
	{
	double sigSquare_i = 0, m_i = 0;
	for ( size_t j = 0; j < NumBinsY; ++j ) 
	  {
	  sigSquare_i += (MathUtil::Square(j) * this->JointBins[ i + j * RealNumBinsX ]);
	  m_i += (j * this->JointBins[ i + j * RealNumBinsX ]);
	  }
	m_i *= (1.0 / ProjectToX(i) );
	(sigSquare_i *= (1.0 / ProjectToX(i))) -= MathUtil::Square(m_i);
	eta += (sigSquare_i * ProjectToX(i));
	}
      }
    
    return eta / (sigSquare * sampleCount); 
  }
  
  /// Copy operator.
  Self& operator=( const Self& other )
  {
    this->BinWidthX = other.BinWidthX;
    this->BinWidthY = other.BinWidthY;
    this->BinOffsetX = other.BinOffsetX;
    this->BinOffsetY = other.BinOffsetY;
    
    this->SetNumBins( other.NumBinsX, other.NumBinsY );
    if ( other.JointBins )
      memcpy( this->JointBins, other.JointBins, RealNumBinsX * RealNumBinsY * sizeof( T ) );
    else
      this->JointBins = NULL;

    return *this;
  }

private:
  /// Real number of X data bins.
  size_t RealNumBinsX;
  
  /// Real number of Y data bins.
  size_t RealNumBinsY;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkJointHistogram_h_included_
