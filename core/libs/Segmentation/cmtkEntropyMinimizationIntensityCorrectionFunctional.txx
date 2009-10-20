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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkThreads.h>
#include <cmtkMemory.h>

#pragma GCC diagnostic ignored "-Wtype-limits"
template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
::UpdateCorrectionFactors()
{
  const int* dims = this->m_InputImage->GetDims();

  // All equation numbers refer to paper by Likar et al., IEEE-TMI 20(12):1398--1410, 2001.
  for ( unsigned int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
    {
    this->m_AddCorrectionAdd[i] = 0;
    this->m_MulCorrectionAdd[i] = 0;
    }

  for ( unsigned int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i )
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
}

#pragma GCC diagnostic ignored "-Wtype-limits"
template<unsigned int NOrderAdd,unsigned int NOrderMul>
typename cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>::ReturnType
cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
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
      if ( dim < PolynomialTypeAdd::NumberOfMonomials )
	this->UpdateBiasFieldAdd();
      else
	this->UpdateBiasFieldMul();
      this->UpdateOutputImage();
      const typename Self::ReturnType upper = this->Evaluate();
      
      v[dim] = v0 - stepScale;
      this->SetParamVector( v );
      if ( dim < PolynomialTypeAdd::NumberOfMonomials )
	this->UpdateBiasFieldAdd();
      else
	this->UpdateBiasFieldMul();
      this->UpdateOutputImage();
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
cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
::UpdateBiasFields( bool foregroundOnly )
{ 
  const size_t numberOfTasks = 2 * ThreadPool::GlobalThreadPool.GetNumberOfThreads();

  ThreadParameters<Self>* taskParameters = Memory::AllocateArray< ThreadParameters<Self> >( numberOfTasks );
  for ( size_t task = 0; task < numberOfTasks; ++task )
    {
    taskParameters[task].thisObject = this;
    }
 
  if ( foregroundOnly )
    ThreadPool::GlobalThreadPool.Run( UpdateBiasFieldsThreadFunc, numberOfTasks, taskParameters );
  else
    ThreadPool::GlobalThreadPool.Run( UpdateBiasFieldsAllThreadFunc, numberOfTasks, taskParameters );
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
::UpdateBiasFieldsThreadFunc( void *args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  ThreadParameters<Self>* threadParameters = static_cast<ThreadParameters<Self>*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;

  const int* dims = ThisConst->m_InputImage->GetDims();
  const UniformVolume* inputImage = ThisConst->m_InputImage;
  TypedArray* outputData = const_cast<TypedArray*>( This->m_OutputImage->GetData().GetPtr() );

  float* biasFieldPtrAdd = This->m_BiasFieldAdd->GetDataPtrTemplate();
  float* biasFieldPtrMul = This->m_BiasFieldMul->GetDataPtrTemplate();

  Types::Coordinate* monomialsVec = This->m_MonomialsVec + (threadIdx * ThisConst->m_MonomialsPerThread);

  const int zFrom = taskIdx * (dims[2] / taskCnt);
  const int zTo = std::max<int>( (taskIdx+1) * (dims[2] / taskCnt), dims[2] );
  
  Types::DataItem value;
  size_t ofs = zFrom * dims[0] * dims[1];
  for ( int z = zFrom; z < zTo; ++z )
    {
    const Types::Coordinate Z = 2.0*(z-dims[2]/2) / dims[2];

    for ( int y = 0; y < dims[1]; ++y )
      {
      const Types::Coordinate Y = 2.0*(y-dims[1]/2) / dims[1];

      for ( int x = 0; x < dims[0]; ++x, ++ofs )
	{
	const Types::Coordinate X = 2.0*(x-dims[0]/2) / dims[0];

	Types::Coordinate normMul = 1.0;
	Types::Coordinate normAdd = 0.0;
	if ( This->m_ForegroundMask[ofs] )
	  {
	  if ( inputImage->GetDataAt( value, ofs ) )
	    {
	    // Normalization according to Eq (13)
	    PolynomialTypeMul::EvaluateAllMonomials( monomialsVec, X, Y, Z );
	    for ( unsigned int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i )
	      {
	      normMul += ThisConst->m_CoefficientsMul[i] * ( monomialsVec[i] - ThisConst->m_AddCorrectionMul[i] );
	      }
	    PolynomialTypeAdd::EvaluateAllMonomials( monomialsVec, X, Y, Z );
	    for ( unsigned int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
	      {
	      normAdd += ThisConst->m_CoefficientsAdd[i] * ( monomialsVec[i] - ThisConst->m_AddCorrectionAdd[i] );
	      }
	    }
	  }
	biasFieldPtrAdd[ofs] = static_cast<float>( normAdd );
	biasFieldPtrMul[ofs] = static_cast<float>( normMul );
	}
      }
    }
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
::UpdateBiasFieldsAllThreadFunc( void *args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  ThreadParameters<Self>* threadParameters = static_cast<ThreadParameters<Self>*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;

  const int* dims = ThisConst->m_InputImage->GetDims();
  const UniformVolume* inputImage = ThisConst->m_InputImage;
  TypedArray* outputData = const_cast<TypedArray*>( This->m_OutputImage->GetData().GetPtr() );

  float* biasFieldPtrAdd = This->m_BiasFieldAdd->GetDataPtrTemplate();
  float* biasFieldPtrMul = This->m_BiasFieldMul->GetDataPtrTemplate();

  Types::Coordinate* monomialsVec = This->m_MonomialsVec + (threadIdx * ThisConst->m_MonomialsPerThread);

  const int zFrom = taskIdx * (dims[2] / taskCnt);
  const int zTo = std::max<int>( (taskIdx+1) * (dims[2] / taskCnt), dims[2] );

  Types::DataItem value;
  size_t ofs = zFrom * dims[0] * dims[1];
  for ( int z = zFrom; z < zTo; ++z )
    {
    const Types::Coordinate Z = 2.0*(z-dims[2]/2) / dims[2];

    for ( int y = 0; y < dims[1]; ++y )
      {
      const Types::Coordinate Y = 2.0*(y-dims[1]/2) / dims[1];

      for ( int x = 0; x < dims[0]; ++x, ++ofs )
	{
	const Types::Coordinate X = 2.0*(x-dims[0]/2) / dims[0];

	Types::Coordinate normMul = 1.0;
	Types::Coordinate normAdd = 0.0;

	if ( inputImage->GetDataAt( value, ofs ) )
	  {
	  // Normalization according to Eq (13)
	  PolynomialTypeMul::EvaluateAllMonomials( monomialsVec, X, Y, Z );
	  for ( unsigned int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i )
	    {
	    normMul += ThisConst->m_CoefficientsMul[i] * ( monomialsVec[i] - ThisConst->m_AddCorrectionMul[i] );
	    }
	  PolynomialTypeAdd::EvaluateAllMonomials( monomialsVec, X, Y, Z );
	  for ( unsigned int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
	    {
	    normAdd += ThisConst->m_CoefficientsAdd[i] * ( monomialsVec[i] - ThisConst->m_AddCorrectionAdd[i] );
	    }
	  }
	biasFieldPtrAdd[ofs] = static_cast<float>( normAdd );
	biasFieldPtrMul[ofs] = static_cast<float>( normMul );
	}
      }
    }
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
::UpdateBiasFieldAdd( const bool foregroundOnly )
{
  const size_t numberOfTasks = 2 * ThreadPool::GlobalThreadPool.GetNumberOfThreads();

  ThreadParameters<Self>* taskParameters = Memory::AllocateArray< ThreadParameters<Self> >( numberOfTasks );
  for ( size_t task = 0; task < numberOfTasks; ++task )
    {
    taskParameters[task].thisObject = this;
    }

  if ( foregroundOnly )
    ThreadPool::GlobalThreadPool.Run( UpdateBiasFieldAddThreadFunc, numberOfTasks, taskParameters );
  else
    ThreadPool::GlobalThreadPool.Run( UpdateBiasFieldAddAllThreadFunc, numberOfTasks, taskParameters );
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
::UpdateBiasFieldAddThreadFunc( void *args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  ThreadParameters<Self>* threadParameters = static_cast<ThreadParameters<Self>*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;

  const int* dims = ThisConst->m_InputImage->GetDims();
  const UniformVolume* inputImage = ThisConst->m_InputImage;
  TypedArray* outputData = const_cast<TypedArray*>( This->m_OutputImage->GetData().GetPtr() );

  float* biasFieldPtrAdd = This->m_BiasFieldAdd->GetDataPtrTemplate();

  Types::Coordinate* monomialsVec = This->m_MonomialsVec + (threadIdx * ThisConst->m_MonomialsPerThread);

  const int zFrom = taskIdx * (dims[2] / taskCnt);
  const int zTo = std::max<int>( (taskIdx+1) * (dims[2] / taskCnt), dims[2] );
  
  size_t ofs = zFrom * dims[0] * dims[1];
  for ( int z = zFrom; z < zTo; ++z )
    {
    const Types::Coordinate Z = 2.0*(z-dims[2]/2) / dims[2];

    for ( int y = 0; y < dims[1]; ++y )
      {
      const Types::Coordinate Y = 2.0*(y-dims[1]/2) / dims[1];

      for ( int x = 0; x < dims[0]; ++x, ++ofs )
	{
	const Types::Coordinate X = 2.0*(x-dims[0]/2) / dims[0];

	Types::Coordinate normAdd = 0.0;
	if ( This->m_ForegroundMask[ofs] )
	  {
	  Types::DataItem value;
	  if ( inputImage->GetDataAt( value, ofs ) )
	    {
	    // Normalization according to Eq (13)
	    PolynomialTypeAdd::EvaluateAllMonomials( monomialsVec, X, Y, Z );
	    for ( unsigned int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
	      {
	      normAdd += ThisConst->m_CoefficientsAdd[i] * ( monomialsVec[i] - ThisConst->m_AddCorrectionAdd[i] );
	      }
	    }
	  }
	biasFieldPtrAdd[ofs] = static_cast<float>( normAdd );
	}
      }
    }
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
::UpdateBiasFieldAddAllThreadFunc( void *args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  ThreadParameters<Self>* threadParameters = static_cast<ThreadParameters<Self>*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;

  const int* dims = ThisConst->m_InputImage->GetDims();
  const UniformVolume* inputImage = ThisConst->m_InputImage;
  TypedArray* outputData = const_cast<TypedArray*>( This->m_OutputImage->GetData().GetPtr() );

  float* biasFieldPtrAdd = This->m_BiasFieldAdd->GetDataPtrTemplate();

  Types::Coordinate* monomialsVec = This->m_MonomialsVec + (threadIdx * ThisConst->m_MonomialsPerThread);

  const int zFrom = taskIdx * (dims[2] / taskCnt);
  const int zTo = std::max<int>( (taskIdx+1) * (dims[2] / taskCnt), dims[2] );

  size_t ofs = zFrom * dims[0] * dims[1];
  for ( int z = zFrom; z < zTo; ++z )
    {
    const Types::Coordinate Z = 2.0*(z-dims[2]/2) / dims[2];

    for ( int y = 0; y < dims[1]; ++y )
      {
      const Types::Coordinate Y = 2.0*(y-dims[1]/2) / dims[1];

      for ( int x = 0; x < dims[0]; ++x, ++ofs )
	{
	const Types::Coordinate X = 2.0*(x-dims[0]/2) / dims[0];

	Types::Coordinate normAdd = 0.0;
	Types::DataItem value;
	if ( inputImage->GetDataAt( value, ofs ) )
	  {
	  // Normalization according to Eq (13)
	  PolynomialTypeAdd::EvaluateAllMonomials( monomialsVec, X, Y, Z );
	  for ( unsigned int i = 0; i < PolynomialTypeAdd::NumberOfMonomials; ++i )
	    {
	    normAdd += ThisConst->m_CoefficientsAdd[i] * ( monomialsVec[i] - ThisConst->m_AddCorrectionAdd[i] );
	    }
	  }
	biasFieldPtrAdd[ofs] = static_cast<float>( normAdd );
	}
      }
    }
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
::UpdateBiasFieldMul( const bool foregroundOnly )
{
  const size_t numberOfTasks = 2 * ThreadPool::GlobalThreadPool.GetNumberOfThreads();

  ThreadParameters<Self>* taskParameters = Memory::AllocateArray< ThreadParameters<Self> >( numberOfTasks );
  for ( size_t task = 0; task < numberOfTasks; ++task )
    {
    taskParameters[task].thisObject = this;
    }

  if ( foregroundOnly )
    ThreadPool::GlobalThreadPool.Run( UpdateBiasFieldMulThreadFunc, numberOfTasks, taskParameters );
  else
    ThreadPool::GlobalThreadPool.Run( UpdateBiasFieldMulAllThreadFunc, numberOfTasks, taskParameters );
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
::UpdateBiasFieldMulThreadFunc( void *args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  ThreadParameters<Self>* threadParameters = static_cast<ThreadParameters<Self>*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;

  const int* dims = ThisConst->m_InputImage->GetDims();
  const UniformVolume* inputImage = ThisConst->m_InputImage;
  TypedArray* outputData = const_cast<TypedArray*>( This->m_OutputImage->GetData().GetPtr() );

  float* biasFieldPtrMul = This->m_BiasFieldMul->GetDataPtrTemplate();

  Types::Coordinate* monomialsVec = This->m_MonomialsVec + (threadIdx * ThisConst->m_MonomialsPerThread);

  const int zFrom = taskIdx * (dims[2] / taskCnt);
  const int zTo = std::max<int>( (taskIdx+1) * (dims[2] / taskCnt), dims[2] );

  size_t ofs = zFrom * dims[0] * dims[1];
  for ( int z = zFrom; z < zTo; ++z )
    {
    const Types::Coordinate Z = 2.0*(z-dims[2]/2) / dims[2];

    for ( int y = 0; y < dims[1]; ++y )
      {
      const Types::Coordinate Y = 2.0*(y-dims[1]/2) / dims[1];

      for ( int x = 0; x < dims[0]; ++x, ++ofs )
	{
	const Types::Coordinate X = 2.0*(x-dims[0]/2) / dims[0];

	Types::Coordinate normMul = 1.0;
	if ( This->m_ForegroundMask[ofs] )
	  {
	  Types::DataItem value;
	  if ( inputImage->GetDataAt( value, ofs ) )
	    {
	    // Normalization according to Eq (13)
	    PolynomialTypeMul::EvaluateAllMonomials( monomialsVec, X, Y, Z );
	    for ( unsigned int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i )
	      {
	      normMul += ThisConst->m_CoefficientsMul[i] * ( monomialsVec[i] - ThisConst->m_AddCorrectionMul[i] );
	      }
	    }
	  }
	biasFieldPtrMul[ofs] = static_cast<float>( normMul );
	}
      }
    }
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctional<NOrderAdd,NOrderMul>
::UpdateBiasFieldMulAllThreadFunc( void *args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  ThreadParameters<Self>* threadParameters = static_cast<ThreadParameters<Self>*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;

  const int* dims = ThisConst->m_InputImage->GetDims();
  const UniformVolume* inputImage = ThisConst->m_InputImage;
  TypedArray* outputData = const_cast<TypedArray*>( This->m_OutputImage->GetData().GetPtr() );

  float* biasFieldPtrMul = This->m_BiasFieldMul->GetDataPtrTemplate();

  Types::Coordinate* monomialsVec = This->m_MonomialsVec + (threadIdx * ThisConst->m_MonomialsPerThread);

  const int zFrom = taskIdx * (dims[2] / taskCnt);
  const int zTo = std::max<int>( (taskIdx+1) * (dims[2] / taskCnt), dims[2] );

  size_t ofs = zFrom * dims[0] * dims[1];
  for ( int z = zFrom; z < zTo; ++z )
    {
    const Types::Coordinate Z = 2.0*(z-dims[2]/2) / dims[2];

    for ( int y = 0; y < dims[1]; ++y )
      {
      const Types::Coordinate Y = 2.0*(y-dims[1]/2) / dims[1];

      for ( int x = 0; x < dims[0]; ++x, ++ofs )
	{
	const Types::Coordinate X = 2.0*(x-dims[0]/2) / dims[0];

	Types::Coordinate normMul = 1.0;
	Types::DataItem value;
	if ( inputImage->GetDataAt( value, ofs ) )
	  {
	  // Normalization according to Eq (13)
	  PolynomialTypeMul::EvaluateAllMonomials( monomialsVec, X, Y, Z );
	  for ( unsigned int i = 0; i < PolynomialTypeMul::NumberOfMonomials; ++i )
	    {
	    normMul += ThisConst->m_CoefficientsMul[i] * ( monomialsVec[i] - ThisConst->m_AddCorrectionMul[i] );
	    }
	  }
	biasFieldPtrMul[ofs] = static_cast<float>( normMul );
	}
      }
    }
}

