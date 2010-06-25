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

#include <cmtkPolynomial.h>

size_t
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice
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
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice
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
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice
::SetInputImage( UniformVolume::SmartConstPtr& inputImage )
{
  this->Superclass::SetInputImage( inputImage );
  this->m_InputImageDevice = UniformVolumeDevice::Create( *inputImage, 512 );
  this->m_NumberOfPixels = inputImage->GetNumberOfPixels();
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice
::SetForegroundMask( const UniformVolume& foregroundMask )
{
  this->Superclass::SetForegroundMask( foregroundMask );
  this->m_ForegroundMaskDevice = UniformVolumeDevice::Create( foregroundMask, 512 );
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice
::SetBiasFieldAdd( const UniformVolume& biasFieldAdd )
{
  this->Superclass::SetBiasFieldAdd( biasFieldAdd );
  this->m_BiasFieldAddDevice->CopyToDevice( this->m_BiasFieldAdd->GetDataPtrTemplate(), this->m_BiasFieldAdd->GetDataSize() );
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice
::SetBiasFieldMul( const UniformVolume& biasFieldMul )
{
  this->Superclass::SetBiasFieldMul( biasFieldMul );
  this->m_BiasFieldMulDevice->CopyToDevice( this->m_BiasFieldMul->GetDataPtrTemplate(), this->m_BiasFieldMul->GetDataSize() );
}
  
void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice
::UpdateCorrectionFactors()
{
#if 0
  const DataGrid::IndexType& dims = this->m_InputImage->GetDims();

  // All equation numbers refer to paper by Likar et al., IEEE-TMI 20(12):1398--1410, 2001.
  for ( unsigned int i = 0; i < this->GetNumberOfMonomialsAdd(); ++i )
    {
    this->m_AddCorrectionAdd[i] = 0;
    this->m_MulCorrectionAdd[i] = 0;
    }

  for ( unsigned int i = 0; i < this->GetNumberOfMonomialsMul(); ++i )
    {
    this->m_AddCorrectionMul[i] = 0;
    this->m_MulCorrectionMul[i] = 0;
    }

  double totalImageEnergy = 0.0;
  size_t foregroundNumberOfPixels = 0;

  // first, compute additive correction factors according to
  // Eqs. (A8) and (A12).
  size_t ofs = 0;
  for ( int z = 0; z < dims[2]; ++z )
    {
    const Types::Coordinate Z = 2.0*(z-dims[2]/2) / dims[2];

    for ( int y = 0; y < dims[1]; ++y )
      {
      const Types::Coordinate Y = 2.0*(y-dims[1]/2) / dims[1];

      for ( int x = 0; x < dims[0]; ++x, ++ofs )
	{
	const Types::Coordinate X = 2.0*(x-dims[0]/2) / dims[0];

	if ( this->m_ForegroundMask[ofs] )
	  {
	  ++foregroundNumberOfPixels;
	  Types::DataItem value;
	  if ( this->m_InputImage->GetDataAt( value, x, y, z ) )
	    totalImageEnergy += value;
	  else
	    value = 0.0;

	  // Eq. (A8)
	  PolynomialTypeAdd::EvaluateAllMonomials( this->m_MonomialsVec, X, Y, Z );
	  for ( unsigned int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
	    {
	    this->m_AddCorrectionAdd[i] += this->m_MonomialsVec[i];
	    }

	  // Eq. (A12)
	  PolynomialTypeMul::EvaluateAllMonomials( this->m_MonomialsVec, X, Y, Z );
	  for ( unsigned int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i )
	    {
	    this->m_AddCorrectionMul[i] += value * this->m_MonomialsVec[i];
	    }
	  }
	}
      }
    }

  // Normalization according to (A8)
  for ( unsigned int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
    {
    this->m_AddCorrectionAdd[i] /= foregroundNumberOfPixels;
    }
  // Normalization according to (A12)
  for ( unsigned int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i )
    {
    this->m_AddCorrectionMul[i] /= totalImageEnergy;
    }

  // Now, compute multiplicative correction factors according to
  // Eqs. (A14) and (A16).
  ofs = 0;
  for ( int z = 0; z < dims[2]; ++z )
    {
    const Types::Coordinate Z = 2.0*(z-dims[2]/2) / dims[2];

    for ( int y = 0; y < dims[1]; ++y )
      {
      const Types::Coordinate Y = 2.0*(y-dims[1]/2) / dims[1];

      for ( int x = 0; x < dims[0]; ++x, ++ofs )
	{
	const Types::Coordinate X = 2.0*(x-dims[0]/2) / dims[0];

	if ( this->m_ForegroundMask[ofs] )
	  {
	  Types::DataItem value;
	  if ( !this->m_InputImage->GetDataAt( value, x, y, z ) )
	    value = 0.0;

	  // Eq. (A8)
	  PolynomialTypeAdd::EvaluateAllMonomials( this->m_MonomialsVec, X, Y, Z );
	  for ( unsigned int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
	    {
	    this->m_MulCorrectionAdd[i] += fabs( this->m_MonomialsVec[i] - this->m_AddCorrectionAdd[i] );
	    }

	  // Eq. (A12)
	  PolynomialTypeMul::EvaluateAllMonomials( this->m_MonomialsVec, X, Y, Z );
	  for ( unsigned int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i )
	    {
	    this->m_MulCorrectionMul[i] += value * fabs( this->m_MonomialsVec[i] - this->m_AddCorrectionMul[i] );
	    }
	  }
	}
      }
    }

  // Normalization according to (A14)
  for ( unsigned int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
    {
    // invert for speedup of application
    this->m_MulCorrectionAdd[i] = foregroundNumberOfPixels / this->m_MulCorrectionAdd[i];
    this->m_StepSizeAdd[i] = 0.0;
    }
  // Normalization according to (A16)
  for ( unsigned int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i )
    {
    // invert for speedup of application
    this->m_MulCorrectionMul[i] = foregroundNumberOfPixels / this->m_MulCorrectionMul[i];
    this->m_StepSizeMul[i] = 0.0;
    }

  // Finally, compute step scale factors according to Eq. (11).
  ofs = 0;
  for ( int z = 0; z < dims[2]; ++z )
    {
    const Types::Coordinate Z = 2.0*(z-dims[2]/2) / dims[2];

    for ( int y = 0; y < dims[1]; ++y )
      {
      const Types::Coordinate Y = 2.0*(y-dims[1]/2) / dims[1];

      for ( int x = 0; x < dims[0]; ++x, ++ofs )
	{
	const Types::Coordinate X = 2.0*(x-dims[0]/2) / dims[0];

	if ( this->m_ForegroundMask[ofs] )
	  {
	  Types::DataItem value;
	  if ( !this->m_InputImage->GetDataAt( value, x, y, z ) )
	    value = 0.0;

	  // Eq. (A8)
	  PolynomialTypeAdd::EvaluateAllMonomials( this->m_MonomialsVec, X, Y, Z );
	  for ( unsigned int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
	    {
	    this->m_StepSizeAdd[i] += fabs( this->m_MulCorrectionAdd[i] * ( this->m_MonomialsVec[i] - this->m_AddCorrectionAdd[i] ) );
	    }

	  // Eq. (A12)
	  PolynomialTypeMul::EvaluateAllMonomials( this->m_MonomialsVec, X, Y, Z );
	  for ( unsigned int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i )
	    {
	    this->m_StepSizeMul[i] += fabs( value * this->m_MulCorrectionMul[i] * ( this->m_MonomialsVec[i] - this->m_AddCorrectionMul[i] ) );
	    }
	  }
	}
      }
    }
  
  // Normalization according to (11)
  for ( unsigned int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
    {
    // invert for speedup of application
    this->m_StepSizeAdd[i] = foregroundNumberOfPixels / this->m_StepSizeAdd[i];
    }
  // Normalization according to (11)
  for ( unsigned int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i )
    {
    // invert for speedup of application
    this->m_StepSizeMul[i] = foregroundNumberOfPixels / this->m_StepSizeMul[i];
    }
#endif
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice
::UpdateBiasFields( const bool foregroundOnly )
{
  this->UpdateBiasFieldAdd( foregroundOnly );
  this->UpdateBiasFieldMul( foregroundOnly );
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice
::UpdateBiasFieldAdd( const bool foregroundOnly )
{
  if ( !this->m_BiasFieldAddDevice )
    this->m_BiasFieldAddDevice = DeviceMemory<float>::Create( this->m_NumberOfPixels, 512 );
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice
::UpdateBiasFieldMul( const bool foregroundOnly )
{
  if ( !this->m_BiasFieldMulDevice )
    this->m_BiasFieldMulDevice = DeviceMemory<float>::Create( this->m_NumberOfPixels, 512 );
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice
::UpdateOutputImage( const bool foregroundOnly )
{
  if ( !this->m_OutputDataDevice )
    this->m_OutputDataDevice = DeviceMemory<float>::Create( this->m_NumberOfPixels, 512 );

  float* input = this->m_InputImageDevice->GetDataOnDevice().Ptr();
  float* output = this->m_OutputDataDevice->Ptr();
  float* biasAdd = this->m_BiasFieldAddDevice ? this->m_BiasFieldAddDevice->Ptr() : NULL;
  float* biasMul = this->m_BiasFieldAddDevice ? this->m_BiasFieldAddDevice->Ptr() : NULL;

  cmtkEntropyMinimizationIntensityCorrectionFunctionalDeviceUpdateOutputImage( input, output, biasAdd, biasMul, this->m_NumberOfPixels );
}
