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

#ifndef __cmtkHistogram_h_included_
#define __cmtkHistogram_h_included_

#include <cmtkconfig.h>

#include <cmtkHistogramBase.h>
#include <cmtkSmartPtr.h>
#include <cmtkMemory.h>

#include <assert.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Histogram of a distribution with bins of arbitrary types.
 * This template is the base class for one-dimensional histograms that can
 * hold integer or real-valued bins, depending on the template parameter.
 *@param T Template parameter: the type of the histogram bins. Can be integral,
 * or double in case of fractional bins.
 *@param FF Template parameter: boolean double flag. Instantiations need to make
 * sure that this parameter is true whenever the instatiated class implements
 * fractional bins, ie., when T is either double or double, and false otherwise.
 */
template<class T>
class Histogram : 
  /// Inherit some non-template functions.
  public HistogramBase 
{
protected:
  /// Array bins.
  T* Bins;
  
public:
  /// This class.
  typedef Histogram<T> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Bin type.
  typedef T BinType;

  /** Constructor.
   */
  Histogram ( const size_t numBins = 0, const bool reset = true ) 
  {
    this->m_NumBins = numBins;
    if ( this->m_NumBins )
      {
      this->Bins = Memory::AllocateArray<T>( this->m_NumBins );
      if ( reset ) this->Reset();
      }
    else
      this->Bins = NULL;
  }

  /// Copy constructor.
  Histogram ( const Self* other, const bool copyData = true ) 
  {
    this->m_NumBins = other->m_NumBins;
    this->m_BinWidth = other->m_BinWidth;
    this->Bins = Memory::AllocateArray<T>( this->m_NumBins );
    
    if ( copyData )
      memcpy( Bins, other->Bins, m_NumBins * sizeof( T ) );
    else
      this->Reset();
  }

  /** Destructor.
   * All bin arrays and the precomputed data bin index arrays are
   * de-allocated.
   */
  virtual ~Histogram () 
  { 
    if (this->Bins) delete[] this->Bins; 
  }

  /// Resize and allocate histogram bins.
  void Resize( const size_t numberOfBins, const bool reset = true )
  {
    if ( numberOfBins != this->m_NumBins )
      {
      if ( this->Bins ) 
	delete[] this->Bins;      

      this->m_BinWidth = this->m_BinWidth * this->m_NumBins / numberOfBins;
      this->m_NumBins = numberOfBins;
      Bins = Memory::AllocateArray<T>( this->m_NumBins );
      }

    if ( reset ) 
      this->Reset();
  }

  /// Copy another histogram without range checking.
  void CopyUnsafe ( const Self& other ) 
  {
    memcpy( Bins, other.Bins, m_NumBins * sizeof( T ) );
  }

  /** Make an identical copy of this object.
   *@param copyData If this is non-zero (default), then the current values in
   * the histogram bins are also duplicated. If zero, only the histogram
   * structure is duplicated.
   *@return The newly created copy of this object.
   */
  Self *Clone ( const bool copyData = true ) const;

  /** Reset computation.
   * This function has to be called before any other computation made with an
   * object of this class. All bin counters are reset to zero, therefore
   * Reset() must also be called before any new computation performed using an
   * already previously used object.
   */
  void Reset ()
  {
    memset( Bins, 0, m_NumBins * sizeof( T ) );
  }

  /// Return number of values in a given bin.
  T GetBin ( const size_t index ) const 
  {
    assert( index < this->m_NumBins );
    return Bins[index];
  }

  /// Set value of a bin.
  void SetBin ( const size_t index, const T value ) 
  {
    assert( index < this->m_NumBins );
    Bins[index] = value;
  }

  /// Return number of values in a given bin using [] operator.
  const T operator[] ( const size_t index ) const 
  {
    assert( index < this->m_NumBins );
    return Bins[index];
  }

  /// Return reference to given bin.
  T& operator[] ( const size_t index ) 
  {
    assert( index < this->m_NumBins );
    return Bins[index];
  }

  /// Set histogram values from existing data array.
  void SetBins( const T* bins ) 
  {
    memcpy( Bins, bins, m_NumBins * sizeof( *Bins ) );
  }

  /// Get histogram values and put them into data array.
  T* GetBins() const 
  {
    T *bins = Memory::AllocateArray<T>( m_NumBins );
    memcpy( bins, Bins, m_NumBins * sizeof( *bins ) );
    return bins;
  }

  /// Get constant pointer to internal histogram values array.
  const T* GetRawBins() const 
  {
    return Bins;
  }

  /** Return total number of samples stored in the histogram.
   */
  T SampleCount () const 
  {
    T sampleCount = 0;
    
    for ( size_t i=0; i<m_NumBins; ++i )
      sampleCount += Bins[i];
    
    return sampleCount;
  }

  /** Return index of bin with highest value.
   */
  size_t GetMaximumBinIndex () const;

  /** Return maximum number of samples stored in any bin.
   */
  T GetMaximumBinValue () const 
  {
    return Bins[ this->GetMaximumBinIndex() ];
  }

  /** Compute entropy of distribution.
   * From the bin counts, the entropy of the distribution of values is 
   * estimated.
   */
  double GetEntropy() const;

  /** Get Kullback-Leibler divergence to other histogram. 
   * The second histogram must have the same number of bins, because the function
   * assumes bin-to-bin correspondence between the two distributions.
   *
   *\note Any bin value ranges set in derived classes are ignored here!
   */
  double GetKullbackLeiblerDivergence( const Self& other ) const;
  
  /** Increment the value of a histogram bin by 1.
   * The histogram field to increment is identified directly by its index;
   * no value-rescaling is done internally.
   *@param sample Index of histogram field.
   */
  void Increment ( const size_t sample ) 
  {
    assert( (sample+1) && (sample<m_NumBins) );

    ++Bins[sample];
  }

  /** Add weighted symmetric kernel to histogram.
   *@param sample Index of histogram field.
   */
  void AddWeightedSymmetricKernel( const size_t bin, const size_t kernelRadius, const T* kernel, const T factor = 1 );
  
  /** Add weighted symmetric kernel to histogram.
   *@param bin Index of histogram field.
   */
  void AddWeightedSymmetricKernelFractional( const double bin, const size_t kernelRadius, const T* kernel, const T factor = 1 );

  /** Increment the value of a histogram bins by fractions of 1.
   * The histogram field to increment is identified directly by its index;
   * no value-rescaling is done internally. The index for this function can be
   * a fractional value, in which case the entry is linearly distributed among
   * neighbouring bins.
   *\note If the bin type of this template object is an integer type, then
   * only the lower of two candidate bins will be decremented by 1.
   *@param sample Fractional index of histogram field.
   */
  void IncrementFractional ( const double bin ) 
  {
    assert( (bin >= 0) && (bin < m_NumBins) );

    const T relative = static_cast<T>( bin - floor(bin) );
    Bins[static_cast<size_t>(bin)] += (1 - relative);
    if ( bin<(m_NumBins-1) )
      Bins[static_cast<size_t>(bin+1)] += relative;
  }

  /** Increment the value of a histogram bin by a given weight.
   * The histogram field to increment is identified directly by its index;
   * no value-rescaling is done internally.
   *@param sample Index of histogram field.
   *@param weight Weight of the current value, i.e., real value that the given
   * bin is incremented by.
   */
  void Increment ( const size_t sample, const double weight ) 
  {
    assert( (sample+1) && (sample<m_NumBins) );

    Bins[sample] += static_cast<T>( weight );
  }

  /** Decrement the value of a histogram bin by 1.
   * The histogram field to decrement is identified directly by its index;
   * no value-rescaling is done internally. Make sure that a value has actually
   * been added to this bin before - otherwise, the next entropy computation my
   * give some unexpected results.
   *@param sample Index of histogram field in direction.
   */
  void Decrement ( const size_t sample ) 
  {
    assert( (sample+1) && (sample<m_NumBins) );
    
    assert( Bins[sample] >= 1 );
    --Bins[sample];
  }

  /** Decrement the value of a histogram bins by fractions of 1.
   * The histogram field to increment is identified directly by its index;
   * no value-rescaling is done internally. The index for this function can be
   * a fractional value, in which case the entry is linearly distributed among
   * neighbouring bins.
   *\note If the bin type of this template object is an integer type, then
   * only the lower of two candidate bins will be decremented by 1.
   *@param sample Fractional index of histogram field.
   */
  void DecrementFractional ( const double bin ) 
  {
    assert( (bin >= 0) && (bin < m_NumBins) );

    T relative = static_cast<T>( bin - floor(bin) );
    Bins[static_cast<size_t>(bin)] -= (1 - relative);
    if ( bin<(m_NumBins-1) )
      Bins[static_cast<size_t>(bin+1)] -= relative;
  }

  /** Decrement the value of a histogram bin by given weight.
   * The histogram field to decrement is identified directly by its index;
   * no value-rescaling is done internally. Make sure that a value has actually
   * been added to this bin before - otherwise, the next entropy computation my
   * give some unexpected results.
   *@param sample Index of histogram field in direction.
   *@param weight Weight of the current value, i.e., real value that the given
   * bin is decremented by.
   */
  void Decrement ( const size_t sample, const double weight ) 
  {
    assert( (sample+1) && (sample<m_NumBins) );
    
    assert( Bins[sample] >= weight );
    Bins[sample] -= static_cast<T>( weight );
  }

  /** Add values from another histogram.
   * Adding is done by corresponding bins. The caller has to make sure that
   * both histograms actually have the same number and arrangement of bins.
   * It is also a good idea to ensure that the data range of these bins is
   * the same in both objects. Both can be guaranteed if one histogram was
   * created from the other by a call to Clone() for example.
   *@param other A pointer to the other histogram. Its bin values are added to
   * this object's bins.
   */
  void AddHistogram ( const Self& other );

  /** Subtract bin values from another histogram.
   * Subtraction is done by corresponding bins. The caller has to make sure 
   * that both histograms actually have the same number of 
   * bins. It is also a good idea to ensure that the data ranges of these bins
   * are the same in both objects. Both can be guaranteed if one histogram was
   * created from the other by a call to Clone() for example.
   *@param other A pointer to the other histogram. Its bin values are
   * subtracted this object's bins.
   */
  void RemoveHistogram ( const Self& other );


  /// Convert this histogram to a cumulative histogram (in place).
  void ConvertToCumulative()
  {
    for ( size_t idx = 1; idx < this->m_NumBins; ++idx )
      {
      this->Bins[idx] += this->Bins[idx-1];
      }
  }
  /** Normalize histogram values by their total sum.
   *@param normalizeTo All histogram bins are scaled by a common factor so that
   * their sum matches the value of this parameter.
   */
  void Normalize( const T normalizeTo = 1 );

  /** Normalize histogram values by their maximum.
   *@param normalizeTo All histogram bins are scaled by a common factor so that
   * their maximum matches the value of this parameter.
   */
  void NormalizeMaximum( const T normalizeTo = 1 );

  /** Compute approximate percentile value from histogram.
   */
  Types::DataItem GetPercentile( const Types::DataItem percentile /**!< The percentile to be computed. Value must be between 0 and 1.*/ ) const;
};

//@}

} // namespace cmtk

#include <cmtkHistogram.txx>

#endif // #ifndef __cmtkHistogram_h_included_
