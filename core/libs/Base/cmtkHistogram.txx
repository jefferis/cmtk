/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class T>
void
Histogram<T>
::AddWeightedSymmetricKernel 
( const size_t bin, const size_t kernelRadius, const T* kernel, const T factor ) 
{
  this->m_Bins[bin] += factor * kernel[0];
  for ( size_t idx = 1; idx < kernelRadius; ++idx )
    {
    const T increment = factor * kernel[idx];
    if ( (bin + idx) < this->GetNumBins() )
      this->m_Bins[bin + idx] += increment;
    if ( bin >= idx )
      this->m_Bins[bin - idx] += increment;
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
  
  if ( (binIdx > 0) && (binIdx+1 < this->GetNumBins()) )
    {
    this->m_Bins[binIdx] += (1-relative) * factor * kernel[0];
    this->m_Bins[binIdx+1] += relative * factor * kernel[0];
    }
  
  for ( size_t idx = 1; idx < kernelRadius; ++idx )
    {
    const T increment = factor * kernel[idx];
    
    const size_t upIdx = binIdx+idx+1;
    if ( upIdx < this->GetNumBins() )
      {
      this->m_Bins[upIdx-1] += (1-relative) * increment;
      this->m_Bins[upIdx  ] += relative * increment;
      }
    
    const int dnIdx = binIdx-idx;
    if ( dnIdx >= 0 )
      {
      this->m_Bins[dnIdx  ] += (1-relative) * increment;
      this->m_Bins[dnIdx+1] += relative * increment;
      }
    }
}

template<class T>
double
Histogram<T>
::GetKullbackLeiblerDivergence( const Self& other ) const 
{
  assert( this->GetNumBins() == other.GetNumBins() );

  const T sampleCount = this->SampleCount();
  const T sampleCountOther = other.SampleCount();

  double dKL = 0;
  for ( size_t i=0; i<this->GetNumBins(); ++i ) 
    {
    if ( this->m_Bins[i] ) 
      {
      const double pX = ((double)this->m_Bins[i]) / sampleCount;
      const double qX = ((double)other.m_Bins[i]) / sampleCountOther;
      dKL += pX*log(pX/qX);
      }
    }
  return dKL;
}

} // namespace cmtk
