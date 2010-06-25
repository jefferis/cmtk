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
      this->SetParamVector( v );
      this->UpdateOutputImageDevice();
      const typename Self::ReturnType upper = this->Evaluate();
      
      v[dim] = v0 - stepScale;
      this->SetParamVector( v );
      this->UpdateOutputImageDevice();
      const  typename Self::ReturnType lower = this->Evaluate();
      
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
      parameters[i] = static_cast<float>( this->m_CoefficientsMul[i] );
      corrections[i] = static_cast<float>( this->m_AddCorrectionMul[i] );
      }
    cmtkEntropyMinimizationIntensityCorrectionFunctionalDeviceUpdateOutputImage( output, input, dims0, dims1, dims2, NOrderAdd, 0 /*multiply*/, Self::PolynomialTypeAdd::NumberOfMonomials, &parameters[0], &corrections[0] );
    }

  this->m_OutputDataDevice->CopyFromDevice( this->m_OutputImage->GetData()->GetDataPtr(), this->m_NumberOfPixels );
}
