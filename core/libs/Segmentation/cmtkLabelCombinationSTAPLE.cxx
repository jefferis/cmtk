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

#include <cmtkLabelCombinationSTAPLE.h>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

LabelCombinationSTAPLE::LabelCombinationSTAPLE( const std::vector<TypedArray::SmartPtr>& data, const int maxIterations, const ScalarDataType resultType )
{
  const size_t numberOfInputs = data.size();
  const size_t numberOfPixels = data[ 0 ]->GetDataSize();
  this->m_Result = TypedArray::SmartPtr( TypedArray::Create( resultType, numberOfPixels ) );

  // compute initial estimate as the average of all inputs;
  // this is also the first E-step with all p/q equal to 0.5
  float totalSum = 0;
#pragma omp parallel for reduction(+:totalSum)
  for ( size_t n = 0; n < numberOfPixels; ++n )
    {
    Types::DataItem w = 0;
    for ( size_t i = 0; i < numberOfInputs; ++i )
      {
      Types::DataItem value;
      if ( data[i]->Get( value, n ) )
	{
	w += value;
	totalSum += value;
	}
      }
    this->m_Result->Set( w / numberOfInputs, n );
    }

  // global prior probability
  const float globalPrior = totalSum / (numberOfInputs * numberOfPixels );

  // expert parameters
  this->m_VecP.resize( numberOfInputs );
  this->m_VecQ.resize( numberOfInputs );

  // iterate
  for ( int it = 0; it < maxIterations; ++it )
    {
    // M-step
    for ( size_t i = 0; i < numberOfInputs; ++i ) 
      {
      this->m_VecP[i] = this->m_VecQ[i] = 0;
      }

    float sumW = 0;
    for ( size_t n = 0; n < numberOfPixels; ++n )
      {
      Types::DataItem w;
      this->m_Result->Get( w, n );
      sumW += w;
      
      for ( size_t i = 0; i < numberOfInputs; ++i ) 
	{
	Types::DataItem value;
	data[i]->Get( value, n );
	this->m_VecP[i] += w * value;
	this->m_VecQ[i] += (1.0 - w) * (1.0 - value);
	}
      }

    for ( size_t i = 0; i < numberOfInputs; ++i ) 
      {
      this->m_VecP[i] /= sumW;
      this->m_VecQ[i] /= (numberOfPixels - sumW);
      }

    // E-step
#pragma omp parallel for
    for ( size_t n = 0; n < numberOfPixels; ++n )
      {
      float alpha = globalPrior;
      float beta = (1.0-globalPrior);
    
      Types::DataItem w = 0;
      for ( size_t i = 0; i < numberOfInputs; ++i )
	{
	if ( data[i]->Get( w, n ) )
	  {
	  alpha *= (1-w-m_VecP[i]);
	  beta *= (w-m_VecQ[i]);
	  }
	}
      this->m_Result->Set( alpha / (alpha+beta), n );
      }
    }
}

} // namespace cmtk
