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

#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice.h"
#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice_kernels.h"

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice<NOrderAdd,NOrderMul>
::SetInputImage( UniformVolume::SmartConstPtr& inputImage )
{
  this->Superclass::SetInputImage( inputImage );
  this->m_InputImageDevice = DeviceUniformVolume::Create( *inputImage, 512 );
  this->m_NumberOfPixels = inputImage->GetNumberOfPixels();

  this->m_HistogramDevice = DeviceHistogram::Create( this->m_NumberOfHistogramBins );
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice<NOrderAdd,NOrderMul>
::SetForegroundMask( const UniformVolume& foregroundMask )
{
  this->Superclass::SetForegroundMask( foregroundMask );

  std::vector<int> maskCopy( this->m_NumberOfPixels );
  for ( size_t i = 0; i < this->m_NumberOfPixels; ++i )
    {
    if ( this->m_ForegroundMask[i] )
      maskCopy[i] = 1;
    else
      maskCopy[i] = 0;
    }

  this->m_ForegroundMaskDevice = DeviceMemory<int>::Create( this->m_NumberOfPixels, &maskCopy[0], 512 );
}

#pragma GCC diagnostic ignored "-Wtype-limits"
template<unsigned int NOrderAdd,unsigned int NOrderMul>
typename cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice<NOrderAdd,NOrderMul>::ReturnType
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice<NOrderAdd,NOrderMul>
::EvaluateWithGradient
( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step )
{ 
  const typename Self::ReturnType baseValue = this->EvaluateAt( v );
  
  for ( size_t dim = 0; dim < this->VariableParamVectorDim(); ++dim ) 
    {
    const Types::Coordinate stepScale = this->GetParamStep( dim, step );
    if ( stepScale <= 0 ) 
      {
      g[dim] = 0;
      } 
    else
      {
      const Types::Coordinate v0 = v[dim];
      
      v[dim] += stepScale;
      const typename Self::ReturnType upper = this->EvaluateAt( v );
      
      v[dim] = v0 - stepScale;
      const  typename Self::ReturnType lower = this->EvaluateAt( v );
      
      v[dim] = v0;
      
      if ( (upper > baseValue) || (lower > baseValue) ) 
	{
	g[dim] = upper-lower;
	} 
      else 
	{
	g[dim] = 0;
	}
      }
    }  

  return baseValue;
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice<NOrderAdd,NOrderMul>
::UpdateOutputImageDevice()
{
  if ( !this->m_OutputDataDevice )
    this->m_OutputDataDevice = DeviceMemory<float>::Create( this->m_NumberOfPixels, 512 );

  float* input = this->m_InputImageDevice->GetDataOnDevice().Ptr();
  float* output = this->m_OutputDataDevice->Ptr();

  const int dims0 = this->m_InputImage->m_Dims[0];
  const int dims1 = this->m_InputImage->m_Dims[1];
  const int dims2 = this->m_InputImage->m_Dims[2];

  if ( Self::PolynomialTypeMul::NumberOfMonomials )
    {
    std::vector<float> parameters( Self::PolynomialTypeMul::NumberOfMonomials ), corrections( Self::PolynomialTypeMul::NumberOfMonomials );
    for ( size_t i = 0; i < Self::PolynomialTypeMul::NumberOfMonomials; ++i )
      {
      parameters[i] = static_cast<float>( this->m_CoefficientsMul[i] );
      corrections[i] = static_cast<float>( this->m_AddCorrectionMul[i] );
      }
    cmtkEntropyMinimizationIntensityCorrectionFunctionalDeviceUpdateOutputImage( output, input, dims0, dims1, dims2, NOrderMul, 1 /*multiply*/, Self::PolynomialTypeMul::NumberOfMonomials, &parameters[0], &corrections[0] );
    input = output; // if additive bias also, apply to output of multiplicative stage
    }

  if ( Self::PolynomialTypeAdd::NumberOfMonomials )
    {
    std::vector<float> parameters( Self::PolynomialTypeAdd::NumberOfMonomials ), corrections( Self::PolynomialTypeAdd::NumberOfMonomials );
    for ( size_t i = 0; i < Self::PolynomialTypeAdd::NumberOfMonomials; ++i )
      {
      parameters[i] = static_cast<float>( this->m_CoefficientsAdd[i] );
      corrections[i] = static_cast<float>( this->m_AddCorrectionAdd[i] );
      }
    cmtkEntropyMinimizationIntensityCorrectionFunctionalDeviceUpdateOutputImage( output, input, dims0, dims1, dims2, NOrderAdd, 0 /*multiply*/, Self::PolynomialTypeAdd::NumberOfMonomials, &parameters[0], &corrections[0] );
    }
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
typename cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice<NOrderAdd,NOrderMul>::ReturnType 
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice<NOrderAdd,NOrderMul>
::EvaluateDevice()
{
  const Types::DataItemRange range = this->m_EntropyHistogram->GetRange();
  this->m_HistogramDevice->Reset();
  this->m_HistogramDevice->Populate( *this->m_OutputDataDevice, *this->m_ForegroundMaskDevice, range.m_LowerBound, range.m_UpperBound );

  const float entropy = this->m_HistogramDevice->GetEntropy();

  std::cerr << entropy << std::endl;

  return -entropy;
}
