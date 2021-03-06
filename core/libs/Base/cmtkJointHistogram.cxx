/*
//
//  Copyright 2016 Google, Inc.
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010, 2013 SRI International
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

#include "cmtkJointHistogram.h"

#include <Base/cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class T> void 
JointHistogram<T>::GetMarginalEntropies ( double& HX, double& HY ) 
  const 
{
  const T sampleCount = this->SampleCount();

  if ( sampleCount > 0 )
    {
    HX = HY = 0;
    for ( size_t i=0; i<NumBinsX; ++i ) 
      {
      const double project = this->ProjectToX( i );
      if ( project ) 
	{
	const double pX = project / sampleCount;
	HX -= pX * log(pX);
	}
      }
  
    for ( size_t j=0; j<NumBinsY; ++j ) 
      {
      const double project = this->ProjectToY( j );
      if ( project ) 
	{
	const double pY = project / sampleCount;
	HY -= pY * log(pY);
	}
      }
    }
  else
    {
    HX = HY = 0.0;
    }
}

template<class T> double
JointHistogram<T>::GetJointEntropy() const 
{
  double HXY = 0;
  
  const T sampleCount = this->SampleCount();
  if ( sampleCount > 0 )
    {
    for ( size_t idx = 0; idx < this->m_TotalNumberOfBins; ++idx )
      {
      if ( JointBins[idx] ) 
	{
	const double pXY = ((double)JointBins[idx]) / sampleCount;
	HXY -= pXY * log(pXY);
	}
      }
    }

  return HXY;
}

template<class T> 
Histogram<T>* 
JointHistogram<T>::GetMarginalX() const 
{
  Histogram<T>* marg = new Histogram<T>( NumBinsX );
  
  marg->SetRange( this->GetRangeX() );
  for ( size_t i = 0; i < NumBinsX; ++i )
    (*marg)[i] = this->ProjectToX( i );
  
  return marg;
}

template<class T> 
Histogram<T>* 
JointHistogram<T>::GetMarginalY() const 
{
  Histogram<T>* marg = new Histogram<T>( NumBinsY );
  
  marg->SetRange( this->GetRangeY() );
  for ( size_t i = 0; i < NumBinsY; ++i )
    (*marg)[i] = this->ProjectToY( i );
  
  return marg;
}

template class JointHistogram<int>;
template class JointHistogram<unsigned int>;
template class JointHistogram<long long int>;
template class JointHistogram<float>;
template class JointHistogram<double>;

} // namespace cmtk
