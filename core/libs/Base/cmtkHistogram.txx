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

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Add weighted symmetric kernel to histogram.
   *@param sample Index of histogram field.
   */
template<class T>
void
Histogram<T>
::AddWeightedSymmetricKernel 
( const size_t bin, const size_t kernelRadius, const T* kernel, const T factor ) 
{
//  assert( (0 <= bin) && (bin < this->m_NumBins) );
  this->Bins[bin] += factor * kernel[0];
  for ( size_t idx = 1; idx < kernelRadius; ++idx )
    {
    const T increment = factor * kernel[idx];
    if ( (bin + idx) < this->m_NumBins )
      this->Bins[bin + idx] += increment;
    if ( bin >= idx )
      this->Bins[bin - idx] += increment;
    }
}

template<class T>
void
Histogram<T>
::AddWeightedSymmetricKernelFractional
( const double bin, const size_t kernelRadius, const T* kernel, const T factor ) 
{
  const T relative = static_cast<T>( bin - floor(bin) );
  const size_t binIdx = static_cast<size_t>( bin );
  
  if ( (binIdx > 0) && (binIdx+1 < this->m_NumBins) )
    {
    this->Bins[binIdx] += (1-relative) * factor * kernel[0];
    this->Bins[binIdx+1] += relative * factor * kernel[0];
    }
  
  for ( size_t idx = 1; idx < kernelRadius; ++idx )
    {
    const T increment = factor * kernel[idx];
    
    const size_t upIdx = binIdx+idx+1;
    if ( upIdx < this->m_NumBins )
      {
      this->Bins[upIdx-1] += (1-relative) * increment;
      this->Bins[upIdx  ] += relative * increment;
      }
    
    const int dnIdx = binIdx-idx;
    if ( dnIdx >= 0 )
      {
      this->Bins[dnIdx  ] += (1-relative) * increment;
      this->Bins[dnIdx+1] += relative * increment;
      }
    }
}

template<class T>
double
Histogram<T>
::GetKullbackLeiblerDivergence( const Self& other ) const 
{
  assert( this->m_NumBins == other.m_NumBins );

  const T sampleCount = this->SampleCount();
  const T sampleCountOther = other.SampleCount();

  double dKL = 0;
  for ( size_t i=0; i<this->m_NumBins; ++i ) 
    {
    if ( this->Bins[i] ) 
      {
      const double pX = ((double)this->Bins[i]) / sampleCount;
      const double qX = ((double)other.Bins[i]) / sampleCountOther;
      dKL += pX*log(pX/qX);
      }
    }
  return dKL;
}

} // namespace cmtk
