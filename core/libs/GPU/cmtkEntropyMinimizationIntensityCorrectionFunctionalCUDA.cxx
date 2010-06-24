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

#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDA.h"
#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDA_functions.h"

#include "cmtkUniformVolumeCUDA.h"
#include <cmtkPolynomial.h>

size_t
cmtk::EntropyMinimizationIntensityCorrectionFunctionalCUDA
::GetNumberOfMonomialsAdd() const 
{
  switch ( this->m_PolyDegreeAdd )
    {
    case 1:
      return Polynomial<1,Types::Coordinate>::NumberOfMonomials;
    case 2:
      return Polynomial<2,Types::Coordinate>::NumberOfMonomials;
    case 3:
      return Polynomial<3,Types::Coordinate>::NumberOfMonomials;
    case 4:
      return Polynomial<4,Types::Coordinate>::NumberOfMonomials;
    default:
      break;
    }
  return 0;
}

size_t
cmtk::EntropyMinimizationIntensityCorrectionFunctionalCUDA
::GetNumberOfMonomialsMul() const 
{
  switch ( this->m_PolyDegreeMul )
    {
    case 1:
      return Polynomial<1,Types::Coordinate>::NumberOfMonomials;
    case 2:
      return Polynomial<2,Types::Coordinate>::NumberOfMonomials;
    case 3:
      return Polynomial<3,Types::Coordinate>::NumberOfMonomials;
    case 4:
      return Polynomial<4,Types::Coordinate>::NumberOfMonomials;
    default:
      break;
    }
  return 0;
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalCUDA
::SetInputImage( UniformVolume::SmartConstPtr& inputImage )
{
  this->Superclass::SetInputImage( inputImage );
  this->m_InputImageCUDA = UniformVolumeCUDA::Create( *inputImage );
  this->m_NumberOfPixels = inputImage->GetNumberOfPixels();
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalCUDA
::SetBiasFieldAdd( const UniformVolume& biasFieldAdd )
{
  this->Superclass::SetBiasFieldAdd( biasFieldAdd );
  this->m_BiasFieldAddCUDA->CopyToDevice( this->m_BiasFieldAdd->GetDataPtrTemplate(), this->m_BiasFieldAdd->GetDataSize() );
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalCUDA
::SetBiasFieldMul( const UniformVolume& biasFieldMul )
{
  this->Superclass::SetBiasFieldMul( biasFieldMul );
  this->m_BiasFieldMulCUDA->CopyToDevice( this->m_BiasFieldMul->GetDataPtrTemplate(), this->m_BiasFieldMul->GetDataSize() );
}
  
void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalCUDA
::UpdateCorrectionFactors()
{
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalCUDA
::UpdateBiasFields( const bool foregroundOnly )
{
  this->UpdateBiasFieldAdd( foregroundOnly );
  this->UpdateBiasFieldMul( foregroundOnly );
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalCUDA
::UpdateBiasFieldAdd( const bool foregroundOnly )
{
  if ( !this->m_BiasFieldAddCUDA )
    this->m_BiasFieldAddCUDA = DeviceMemoryCUDA<float>::Create( this->m_NumberOfPixels );
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalCUDA
::UpdateBiasFieldMul( const bool foregroundOnly )
{
  if ( !this->m_BiasFieldMulCUDA )
    this->m_BiasFieldMulCUDA = DeviceMemoryCUDA<float>::Create( this->m_NumberOfPixels );
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalCUDA
::UpdateOutputImage( const bool foregroundOnly )
{
  if ( !this->m_OutputDataCUDA )
    this->m_OutputDataCUDA = DeviceMemoryCUDA<float>::Create( this->m_NumberOfPixels );

  float* input = this->m_InputImageCUDA->GetDataOnDevice().Ptr();
  float* output = this->m_OutputDataCUDA->Ptr();
  float* biasAdd = this->m_BiasFieldAddCUDA ? this->m_BiasFieldAddCUDA->Ptr() : NULL;
  float* biasMul = this->m_BiasFieldAddCUDA ? this->m_BiasFieldAddCUDA->Ptr() : NULL;

  cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImage( input, output, biasAdd, biasMul, this->m_NumberOfPixels );
}
