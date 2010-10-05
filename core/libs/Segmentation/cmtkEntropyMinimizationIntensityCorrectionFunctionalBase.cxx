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

#include "Segmentation/cmtkEntropyMinimizationIntensityCorrectionFunctionalBase.h"

#include "System/cmtkThreadPool.h"
#include "System/cmtkException.h"

#include <algorithm>

#ifdef CMTK_BUILD_DEMO
#  include <IO/cmtkVolumeIO.h>
#endif // #ifdef CMTK_BUILD_DEMO

namespace cmtk
{

/** \addtogroup Segmentation */
//@{
void
EntropyMinimizationIntensityCorrectionFunctionalBase
::SetInputImage( UniformVolume::SmartConstPtr& inputImage )
{
  this->m_InputImage = inputImage;
  this->m_NumberOfPixels = this->m_InputImage->GetNumberOfPixels();
  
  const Types::DataItemRange range = this->m_InputImage->GetData()->GetRange();
  this->m_InputImageRange = range.Width();
  
  if ( this->m_UseLogIntensities )
    {
    this->m_EntropyHistogram = HistogramType::SmartPtr( new LogHistogramType( this->m_NumberOfHistogramBins ) );
    }
  else
    {
    this->m_EntropyHistogram = HistogramType::SmartPtr( new HistogramType( this->m_NumberOfHistogramBins ) );
    }
  // extend value range to accomodate corrected intensities without overflows
  this->m_EntropyHistogram->SetRange( Types::DataItemRange( range.m_LowerBound - this->m_InputImageRange, range.m_UpperBound + this->m_InputImageRange ) );
  
  if ( this->m_ForegroundMask.size() )
    this->UpdateCorrectionFactors();
  
  this->m_BiasFieldAdd = FloatArray::Create( this->m_NumberOfPixels );
  this->m_BiasFieldAdd->Fill( 0.0 );
  this->m_BiasFieldMul = FloatArray::Create( this->m_NumberOfPixels );
  this->m_BiasFieldAdd->Fill( 1.0 );

  this->m_OutputImage = UniformVolume::SmartPtr( this->m_InputImage->CloneGrid() );
  this->m_OutputImage->CreateDataArray( TYPE_FLOAT );
}

void
EntropyMinimizationIntensityCorrectionFunctionalBase
::SetForegroundMask( const UniformVolume& foregroundMask )
{
  const size_t maskPixels = foregroundMask.GetNumberOfPixels();

  if ( maskPixels != this->m_NumberOfPixels )
    {
    throw Exception( "Number of mask pixels does not match number of input image pixels." );
    }

  this->m_ForegroundMask.resize( maskPixels );
  if ( (this->m_SamplingDensity > 0) && (this->m_SamplingDensity < 1) )
    {
    for ( size_t i = 0; i < maskPixels; ++i )
      {
      this->m_ForegroundMask[i] = (foregroundMask.GetDataAt( i ) > 0) && (MathUtil::UniformRandom() <= this->m_SamplingDensity);
      }
    }
  else
    {
    for ( size_t i = 0; i < maskPixels; ++i )
      {
      this->m_ForegroundMask[i] = (foregroundMask.GetDataAt( i ) > 0);
      }
    }
  
  if ( this->m_InputImage )
    this->UpdateCorrectionFactors();
}

UniformVolume::SmartPtr& 
EntropyMinimizationIntensityCorrectionFunctionalBase
::GetOutputImage( CoordinateVector& v, const bool foregroundOnly )
{
  this->SetParamVector( v );
  this->UpdateBiasFields( foregroundOnly );
  this->UpdateOutputImage( foregroundOnly );
  return this->m_OutputImage;
}

EntropyMinimizationIntensityCorrectionFunctionalBase::ReturnType
EntropyMinimizationIntensityCorrectionFunctionalBase
::EvaluateAt( CoordinateVector& v )
{
  this->SetParamVector( v );
  this->UpdateBiasFields();
  this->UpdateOutputImage();
  
#ifdef CMTK_BUILD_DEMO
  static int it = 0;
  this->UpdateOutputImage( false /*foregroundOnly*/ );
  UniformVolume::SmartConstPtr slice = this->m_OutputImage->ExtractSlice( AXIS_Z, this->m_OutputImage->m_Dims[2] / 2 );
  
  char path[PATH_MAX];
  snprintf( path, PATH_MAX, "mrbias-%03d.nii", it++ );
  VolumeIO::Write( *slice, path );
#endif
  
  return this->Evaluate();
}

void
EntropyMinimizationIntensityCorrectionFunctionalBase
::UpdateOutputImage( const bool foregroundOnly )
{
  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  const size_t numberOfTasks = 4 * threadPool.GetNumberOfThreads() - 3;
  
  std::vector<UpdateOutputImageThreadParameters> taskParameters( numberOfTasks );
  for ( size_t task = 0; task < numberOfTasks; ++task )
    {
    taskParameters[task].thisObject = this;
    taskParameters[task].m_ForegroundOnly = foregroundOnly;
    }
  threadPool.Run( UpdateOutputImageThreadFunc, taskParameters );
}
 
void
EntropyMinimizationIntensityCorrectionFunctionalBase
::UpdateOutputImageThreadFunc( void *args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t )
{
  UpdateOutputImageThreadParameters* threadParameters = static_cast<UpdateOutputImageThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  
  const UniformVolume* inputImage = ThisConst->m_InputImage;
  TypedArray::SmartPtr outputData = This->m_OutputImage->GetData();
  const size_t numberOfPixels = inputImage->GetNumberOfPixels();

  const float* biasFieldPtrAdd = ThisConst->m_BiasFieldAdd->GetDataPtrTemplate();
  const float* biasFieldPtrMul = ThisConst->m_BiasFieldMul->GetDataPtrTemplate();

  Types::DataItem value;
  for ( size_t ofs = taskIdx; ofs < numberOfPixels; ofs += taskCnt )
    {
    if ( !threadParameters->m_ForegroundOnly || ThisConst->m_ForegroundMask[ofs] )
      {
      if ( inputImage->GetDataAt( value, ofs ) )
	{
	outputData->Set( value * biasFieldPtrMul[ofs] + biasFieldPtrAdd[ofs], ofs );
	}
      else
	{
	outputData->SetPaddingAt( ofs );
	}
      }	    
    else
      {
      outputData->SetPaddingAt( ofs );
      }
    }
}

} // namespace cmtk
