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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkLogHistogram_h_included_
#define __cmtkLogHistogram_h_included_

#include <cmtkconfig.h>

#include <cmtkHistogram.h>
#include <cmtkSmartPtr.h>
#include <math.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Histogram of log intensities.
 */
template<class T>
class LogHistogram : 
  /// Inherit from non-log histogram..
  public Histogram<T> 
{
public:
  /// This class.
  typedef LogHistogram<T> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Parent class.
  typedef Histogram<T> Superclass;

  /** Constructor.
   */
  LogHistogram ( const size_t numBins = 0, const bool reset = true ) : Histogram<T>( numBins, reset ) {}

  /// Copy constructor.
  LogHistogram ( const Self* other, const bool copyData = true ) : Histogram<T>( other, copyData ) {}

  /** Destructor.
   * All bin arrays and the precomputed data bin index arrays are
   * de-allocated.
   */
  virtual ~LogHistogram() {}

  /** Return bin corresponding to a certain value of the distribution.
   *@param value A value from the distribution.
   *@return The index of the bin corresponding to the given value.
   */
  virtual size_t ValueToBin ( const Types::DataItem value ) const 
  {
    const size_t binIndex = static_cast<size_t>( (this->m_NumBins-1) * log(1 + value - this->m_BinsLowerBound) / log( this->m_BinsUpperBound ) );
    return std::max<size_t>( 0, std::min( this->m_NumBins-1, binIndex ) );
  }

  /** Return fractional bin corresponding to a value of the distribution.
   *@param value A value from the distribution.
   *@return The index of the fractional bin index corresponding to the given 
   * value. This value is an integer if and only if the given value is 
   * identical to the lower bound of a bin.
   */
  virtual Types::DataItem ValueToBinFractional ( const Types::DataItem value ) const 
  {
    const Types::DataItem binIndex = (this->m_NumBins-1) * log(1 + value - this->m_BinsLowerBound) / log( this->m_BinsUpperBound );
    return std::max<Types::DataItem>( 0, std::min<Types::DataItem>( this->m_NumBins-1, binIndex ) );
  }
  
  /** Get value range of a given bin.
   */
  virtual void GetRangeBin( const size_t bin, Types::DataItem& from, Types::DataItem& to ) const 
  {
    from = exp( bin * log( this->m_BinsUpperBound ) / (this->m_NumBins-1) ) + this->m_BinsLowerBound - 1;
    to = exp( (bin+1) * log( this->m_BinsUpperBound ) / (this->m_NumBins-1) ) + this->m_BinsLowerBound - 1;
  }
  
  /** Return center of values represented by a certain bin.
   *@param bin Index of a bin from the distribution.
   *@return Average of upper and lower margin values of the given bin.
   */
  virtual Types::DataItem BinToValue ( const size_t bin ) const 
  {
    return exp( (0.5+bin) * log( this->m_BinsUpperBound ) / (this->m_NumBins-1) ) + this->m_BinsLowerBound - 1;
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkLogHistogram_h_included_
