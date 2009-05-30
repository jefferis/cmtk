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

#include <cmtkSplineWarpCongealingFunctional.h>

#include <cmtkMathUtil.h>
#include <cmtkMatrix.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

SplineWarpCongealingFunctional
::SplineWarpCongealingFunctional()
  : m_ForceZeroSumNoAffine( false )
{
  this->m_ParametersPerXform = 0;
  this->m_WarpFastMode = true;
  this->m_PartialGradientMode = false;
  this->m_PartialGradientThreshold = static_cast<Types::DataItem>( 0.01 );
  this->m_DeactivateUninformativeMode = false;
  this->m_NumberOfActiveControlPoints = 0;
  this->m_JacobianConstraintWeight = 0.0;
  this->m_BendingEnergyWeight = 0.0;
}

SplineWarpCongealingFunctional
::~SplineWarpCongealingFunctional()
{
}

void
SplineWarpCongealingFunctional
::SetTemplateGrid
( UniformVolume::SmartPtr& templateGrid, 
  const int downsample,
  const bool useTemplateData )
{
  this->Superclass::SetTemplateGrid( templateGrid, downsample, useTemplateData );
  this->m_EntropyByPixel.resize( this->m_TemplateNumberOfPixels );
  
  if ( this->m_XformVector.size() )
    {
    for ( size_t i = 0; i < this->m_XformVector.size(); ++i )
      {
      this->m_XformVector[i]->RegisterVolume( this->m_TemplateGrid );	
      }      
    this->UpdateVolumesOfInfluence();
    }

  // clear thread storage because we need to re-initialize these.
  this->m_StaticThreadStorage.resize(0);
}

void
SplineWarpCongealingFunctional
::InitializeXforms
( const Types::Coordinate gridSpacing, std::vector<AffineXform::SmartPtr> initialAffineXformsVector, const bool exactSpacing )
{
  this->m_InitialAffineXformsVector = initialAffineXformsVector;  

  this->m_XformVector.resize( this->m_ImageVector.size() );
  this->m_InitialRotationsVector.resize( this->m_ImageVector.size() );

  for ( size_t i = 0; i < this->m_ImageVector.size(); ++i )
    {
    SplineWarpXform::SmartPtr xform( new SplineWarpXform( this->m_TemplateGrid->Size, gridSpacing, initialAffineXformsVector[i], exactSpacing ) );
    xform->RegisterVolume( this->m_TemplateGrid );
    this->m_XformVector[i] = xform;

    this->m_InitialRotationsVector[i] = AffineXform::SmartPtr( initialAffineXformsVector[i] );

    // create all-zero parameter vector
    CoordinateVector v( initialAffineXformsVector[i]->ParamVectorDim(), 0.0 );
    // copy rotation angles
    for ( size_t p = 3; p < 6; ++p )
      v[p] = initialAffineXformsVector[i]->GetParameter( p );
    // create rotation-only transformation
    this->m_InitialRotationsVector[i]->SetParamVector( v );
    }

  this->m_ParametersPerXform = this->m_XformVector[0]->VariableParamVectorDim();

  this->UpdateVolumesOfInfluence();

  // clear thread storage because we need to re-initialize these.
  this->m_StaticThreadStorage.resize(0);
}

void
SplineWarpCongealingFunctional
::RefineTransformationGrids()
{
  for ( size_t i = 0; i < this->m_XformVector.size(); ++i )
    {
    this->GetXformByIndex(i)->Refine();
    this->m_XformVector[i]->RegisterVolume( this->m_TemplateGrid );
    }

  this->m_ParametersPerXform = this->m_XformVector[0]->VariableParamVectorDim();
  this->UpdateVolumesOfInfluence();

  // clear thread storage because we need to re-initialize these.
  this->m_StaticThreadStorage.resize(0);
}

void
SplineWarpCongealingFunctional
::UpdateProbabilisticSamples()
{
  this->Superclass::UpdateProbabilisticSamples();
}

void
SplineWarpCongealingFunctional
::UpdateStandardDeviationByPixel()
{
  this->Superclass::UpdateStandardDeviationByPixel();
  this->UpdateActiveControlPoints();
  this->UpdateParamStepArray();
}

bool
SplineWarpCongealingFunctional
::UpdateParamStepArray()
{
  bool changed = false;

  this->m_ParamStepArray.resize( this->ParamVectorDim() );

  if ( this->m_DeactivateUninformativeMode &&
       ( this->m_ActiveControlPointFlags.size() == this->m_ParametersPerXform / 3 ) )
    {
    for ( size_t param = 0; param < this->ParamVectorDim(); ++param ) 
      {
      const Types::Coordinate pOld = this->m_ParamStepArray[param];
      this->m_ParamStepArray[param] = this->GetParamStep( param );
      if ( ! this->m_ActiveControlPointFlags[(param%this->m_ParametersPerXform)/3] )
	{
	this->m_ParamStepArray[param] = 0;
	}
      if ( pOld != this->m_ParamStepArray[param] )
	changed = true;
      }
    }
  else
    {
    for ( size_t param = 0; param < this->ParamVectorDim(); ++param ) 
      {
      const Types::Coordinate pOld = this->m_ParamStepArray[param];
      this->m_ParamStepArray[param] = this->GetParamStep( param );
      if ( pOld != this->m_ParamStepArray[param] )
	changed = true;
      }
    }
  
  return changed;
}

void
SplineWarpCongealingFunctional
::UpdateVolumesOfInfluence()
{
  const Vector3D templateFrom( this->m_TemplateGrid->m_Origin );
  const Vector3D templateTo(  this->m_TemplateGrid->m_Origin + this->m_TemplateGrid->Size );
  
  this->m_VolumeOfInfluenceArray.resize( this->m_ParametersPerXform / 3 );

  this->m_MaximumNumberOfPixelsPerLineVOI = 0;
  this->m_MaximumNumberOfPixelsVOI = 0;
  
  const SplineWarpXform* xform0 = this->GetXformByIndex(0);
  for ( size_t param = 0; param < this->m_ParametersPerXform; param += 3 ) 
    { 
    Vector3D fromVOI, toVOI;
    xform0->GetVolumeOfInfluence( param, templateFrom, templateTo, fromVOI, toVOI );

    Rect3D* voi = &this->m_VolumeOfInfluenceArray[param/3];
    this->m_TemplateGrid->GetGridRange( fromVOI, toVOI, *voi );
    
    this->m_MaximumNumberOfPixelsVOI = std::max<size_t>( voi->Size(), this->m_MaximumNumberOfPixelsVOI );
    this->m_MaximumNumberOfPixelsPerLineVOI = std::max<size_t>( voi->endX-voi->startX, this->m_MaximumNumberOfPixelsPerLineVOI );
    }
}
  
void
SplineWarpCongealingFunctional
::UpdateActiveControlPoints()
{
  if ( this->m_DeactivateUninformativeMode )
    {
    const size_t numberOfControlPoints = this->m_VolumeOfInfluenceArray.size();
    
    if ( numberOfControlPoints )
      {
      this->m_ActiveControlPointFlags.resize( numberOfControlPoints );
      this->m_NumberOfActiveControlPoints = 0;
      
      const Vector3D templateFrom( this->m_TemplateGrid->m_Origin );
      const Vector3D templateTo( this->m_TemplateGrid->m_Origin + this->m_TemplateGrid->Size );
      Vector3D fromVOI, toVOI;
      
      const Rect3D* voi = &this->m_VolumeOfInfluenceArray[0];
      for ( size_t cp = 0; cp < numberOfControlPoints; ++cp, ++voi )
	{
	bool active = false;
	for ( int z = voi->startZ; (z < voi->endZ) && !active; ++z ) 
	  {
	  for ( int y = voi->startY; (y < voi->endY) && !active; ++y )
	    {
	    size_t ofs = this->m_TemplateGrid->GetOffsetFromIndex( voi->startX, y, z );
	    for ( size_t x = voi->startX; (x < voi->endX)  && !active; ++x, ++ofs )
	      {
	      if ( this->m_StandardDeviationByPixel[ofs] > 0 )
		{
		active = true;
		}
	      }
	    }
	  }
	this->m_ActiveControlPointFlags[cp] = active;
	if ( active ) ++this->m_NumberOfActiveControlPoints;
	}

      StdErr << "Enabled " << this->m_NumberOfActiveControlPoints << "/" << this->m_ParametersPerXform / 3 << " control points.\n";
      }
    }
  else
    {
    this->m_NumberOfActiveControlPoints = this->m_VolumeOfInfluenceArray.size();
    }
}

SplineWarpCongealingFunctional::ReturnType
SplineWarpCongealingFunctional
::Evaluate()
{
  if ( this->m_NeedsUpdateStandardDeviationByPixel )
    this->UpdateStandardDeviationByPixel();

  const size_t numberOfPixels = this->m_TemplateNumberOfPixels;
#ifdef CMTK_BUILD_MPI
  const size_t pixelsPerNode = (numberOfPixels+this->m_SizeMPI-1) / this->m_SizeMPI;
  this->m_EntropyByPixel.resize( pixelsPerNode * this->m_SizeMPI );
  this->m_EntropyByPixelMPI.resize( pixelsPerNode );
#else
  this->m_EntropyByPixel.resize( numberOfPixels );
#endif

  double entropy = 0;
  unsigned int count = 0;

  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  this->m_ThreadHistograms.resize( numberOfThreads );

  ThreadParameterArray<Self,EvaluateThreadParameters> params( this, numberOfThreads );  
  params.RunInParallel( &EvaluateThread );
  
  // gather partial entropies from threads
  for ( size_t thread = 0; thread < numberOfThreads; ++thread )
    {
    entropy += params[thread].m_Entropy;
    count += params[thread].m_Count;
    }
  
#ifdef CMTK_BUILD_MPI
  double partialEntropy = entropy;
  MPI::COMM_WORLD.Allreduce( &partialEntropy, &entropy, 1, MPI::DOUBLE, MPI::SUM );
  
  unsigned int partialCount = count;
  MPI::COMM_WORLD.Allreduce( &partialCount, &count, 1, MPI::UNSIGNED, MPI::SUM );

  const size_t totalSize = sizeof( this->m_EntropyByPixelMPI[0] ) * this->m_EntropyByPixelMPI.size();
  MPI::COMM_WORLD.Allgather( &this->m_EntropyByPixelMPI[0], totalSize, MPI::CHAR, &this->m_EntropyByPixel[0], totalSize, MPI::CHAR );  
#endif
  
  if ( count )
    {
    const double result = entropy / count;
    double constraint = 0;
    if ( this->m_JacobianConstraintWeight > 0 )
      {
      for ( size_t i = 0; i < this->m_XformVector.size(); ++i )
	{
	constraint += dynamic_cast<const SplineWarpXform*>( this->m_XformVector[i].GetPtr() )->GetJacobianConstraint();
	}
      }
    return result - this->m_JacobianConstraintWeight * constraint;
    }
  else
    return -FLT_MAX;
}

CMTK_THREAD_RETURN_TYPE
SplineWarpCongealingFunctional
::EvaluateThread
( void *args )
{
  EvaluateThreadParameters* threadParameters = static_cast<EvaluateThreadParameters*>( args );
  
  Self* This = threadParameters->thisObject;
  const Self* ThisConst = threadParameters->thisObject;
  const int threadID = threadParameters->ThisThreadIndex;
  const int numberOfThreads = threadParameters->NumberOfThreads;
  
  HistogramType& histogram = This->m_ThreadHistograms[threadID];
  histogram.Resize( ThisConst->m_HistogramBins + 2 * ThisConst->m_HistogramKernelRadiusMax, false /*reset*/ );

  double totalEntropy = 0;
  size_t count = 0;

  const size_t numberOfPixels = ThisConst->m_TemplateNumberOfPixels;
#ifdef CMTK_BUILD_MPI  
  const size_t pixelsPerNode = (numberOfPixels+ThisConst->m_SizeMPI-1) / ThisConst->m_SizeMPI;
  const size_t pixelsPerThread = 1 + (pixelsPerNode / numberOfThreads);
  const size_t pixelFromNode = ThisConst->m_RankMPI * pixelsPerNode;
  const size_t pixelFrom = pixelFromNode + threadID * pixelsPerThread;
  const size_t pixelTo = std::min( numberOfPixels, std::min( pixelFromNode + pixelsPerNode, pixelFrom + pixelsPerThread ) );
  size_t mpiOfs = threadID * pixelsPerThread;
#else
  const size_t pixelsPerThread = (numberOfPixels / numberOfThreads);
  const size_t pixelFrom = threadID * pixelsPerThread;
  const size_t pixelTo = std::min( numberOfPixels, pixelFrom + pixelsPerThread );
#endif
  
  const size_t imagesFrom = ThisConst->m_ActiveImagesFrom;
  const size_t imagesTo = ThisConst->m_ActiveImagesTo;
  const byte paddingValue = ThisConst->m_PaddingValue;
  
  for ( size_t ofs = pixelFrom; ofs < pixelTo; ++ofs )
    {
    histogram.Reset();
    const size_t kernelIdx = ThisConst->m_StandardDeviationByPixel[ofs];
    const size_t kernelRadius = ThisConst->m_HistogramKernelRadius[kernelIdx];
    const HistogramBinType* kernel = ThisConst->m_HistogramKernel[kernelIdx];

    bool fullCount = true;
    
    if ( ThisConst->m_UseTemplateData )
      {
      const byte templateValue = ThisConst->m_TemplateData[ofs];
      if ( (fullCount = (templateValue != paddingValue)) )
	{
	histogram.AddWeightedSymmetricKernel( templateValue, kernelRadius, kernel );
	}
      }

    for ( size_t idx = imagesFrom; (idx < imagesTo) && fullCount; ++idx )
      {
      const byte value = ThisConst->m_Data[idx][ofs];
      if ( value != paddingValue )
	{
	histogram.AddWeightedSymmetricKernel( value, kernelRadius, kernel );
	}
      else
	{
	fullCount = false;
	}
      }
    
    if ( fullCount )
      {
      const double entropy = histogram.GetEntropy();
#ifdef CMTK_BUILD_MPI
      This->m_EntropyByPixelMPI[mpiOfs++] = entropy;
#else
      This->m_EntropyByPixel[ofs] = entropy;
#endif
      totalEntropy -= entropy;
      ++count;
      }
    else
      {
#ifdef CMTK_BUILD_MPI
      This->m_EntropyByPixelMPI[mpiOfs++] = 0;
#else
      This->m_EntropyByPixel[ofs] = 0;
#endif
      }
    }
  
  threadParameters->m_Entropy = totalEntropy;
  threadParameters->m_Count = count;

  return CMTK_THREAD_RETURN_VALUE;
}

void
SplineWarpCongealingFunctional
::InterpolateImage
( const size_t idx, byte* const destination )
{
  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  ThreadParameterArray<Self,InterpolateImageThreadParameters> 
    params( this, numberOfThreads );

  for ( size_t thread = 0; thread < numberOfThreads; ++thread )
    {
    params[thread].m_Idx = idx;    
    params[thread].m_Destination = destination;    
    }
  
  params.RunInParallel( &InterpolateImageThread );
}

CMTK_THREAD_RETURN_TYPE
SplineWarpCongealingFunctional
::InterpolateImageThread
( void* args )
{
  InterpolateImageThreadParameters* threadParameters = static_cast<InterpolateImageThreadParameters*>( args );
  
  const Self* This = threadParameters->thisObject;
  const int threadID = threadParameters->ThisThreadIndex;
  const int numberOfThreads = threadParameters->NumberOfThreads;
  const size_t idx = threadParameters->m_Idx;
  byte* destination = threadParameters->m_Destination;

  const SplineWarpXform* xform = This->GetXformByIndex(idx);
  const UniformVolume* target = This->m_ImageVector[idx];
  const byte* dataPtr = static_cast<const byte*>( target->GetData()->GetDataPtr() );

  Vector3D v;
  byte value;

  const byte paddingValue = This->m_PaddingValue;
  const byte backgroundValue = This->m_UserBackgroundFlag ? This->m_PrivateUserBackgroundValue : paddingValue;

  const int dimsX = This->m_TemplateGrid->GetDims( AXIS_X );
  const int dimsY = This->m_TemplateGrid->GetDims( AXIS_Y );
  const int dimsZ = This->m_TemplateGrid->GetDims( AXIS_Z );

  const int rowCount = ( dimsY * dimsZ );
  const int rowFrom = ( rowCount / numberOfThreads ) * threadID;
  const int rowTo = ( threadID == (numberOfThreads-1) ) ? rowCount : ( rowCount / numberOfThreads ) * ( threadID + 1 );
  int rowsToDo = rowTo - rowFrom;
  
  int yFrom = rowFrom % dimsY;
  int zFrom = rowFrom / dimsY;
  
  byte *wptr = destination + rowFrom * dimsX;
  for ( int z = zFrom; (z < dimsZ) && rowsToDo; ++z ) 
    {
    for ( int y = yFrom; (y < dimsY) && rowsToDo; yFrom = 0, ++y, --rowsToDo )
      {
      for ( int x = 0; x < dimsX; ++x )
	{
	xform->GetTransformedGridNonVirtual( v, x, y, z );
	if ( target->ProbeData( value, dataPtr, v ) )
	  {
	  *wptr = value;
	  }
	else
	  {
	  *wptr = backgroundValue;
	  }
	
	++wptr;
	}
      }
    }

  return CMTK_THREAD_RETURN_VALUE;
}

void
SplineWarpCongealingFunctional
::ForceZeroSumGradient( CoordinateVector& g ) const
{
  const size_t numberOfXforms = this->m_XformVector.size();

  if ( this->m_ForceZeroSumNoAffine )
    {
    for ( size_t xform = 0; xform < numberOfXforms; ++xform )
      {
      Types::Coordinate* gX = &g[xform*this->m_ParametersPerXform];
      const AffineXform* affineXform = this->m_InitialRotationsVector[xform]->GetInverse();
      if ( affineXform )
	{
#pragma omp parallel for
	for ( size_t param = 0; param < this->m_ParametersPerXform; param += 3 )
	  {
	  affineXform->RotateScaleShear( gX+param, gX+param );
	  }
	}
      }
    }
  
  this->Superclass::ForceZeroSumGradient( g );

  if ( this->m_ForceZeroSumNoAffine )
    {
    for ( size_t xform = 0; xform < numberOfXforms; ++xform )
      {
      Types::Coordinate* gX = &g[xform*this->m_ParametersPerXform];
      const AffineXform* affineXform = this->m_InitialRotationsVector[xform];
      if ( affineXform )
	{
#pragma omp parallel for
	for ( size_t param = 0; param < this->m_ParametersPerXform; param += 3 )
	  {
	  affineXform->RotateScaleShear( gX+param, gX+param );
	  }
	}
      }
    }
}

bool
SplineWarpCongealingFunctional
::Wiggle()
{
  bool wiggle = this->Superclass::Wiggle();

  if ( this->m_PartialGradientMode )
    {
    wiggle = wiggle || this->UpdateParamStepArray();
    }

  return wiggle;
}

//@}

} // namespace cmtk

#ifdef CMTK_BUILD_MPI
#  include "cmtkSplineWarpCongealingFunctionalMPI.txx"
#else
#  include "cmtkSplineWarpCongealingFunctional.txx"
#endif
