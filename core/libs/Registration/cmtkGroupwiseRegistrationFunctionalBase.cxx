/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include <Registration/cmtkGroupwiseRegistrationFunctionalBase.h>

#include <Base/cmtkMathUtil.h>
#include <IO/cmtkVolumeIO.h>
#include <System/cmtkConsole.h>
#include <System/cmtkThreadPool.h>

#include <Base/cmtkAnatomicalOrientation.h>

#include <Base/cmtkInterpolator.h>
#include <Base/cmtkTypedArrayFunctionHistogramMatching.h>
#include <Base/cmtkUniformVolumeGaussianFilter.h>
#include <Registration/cmtkReformatVolume.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

GroupwiseRegistrationFunctionalBase
::GroupwiseRegistrationFunctionalBase() 
  : m_FreeAndRereadImages( false ),
    m_ForceZeroSum( false ),
    m_ForceZeroSumFirstN( 0 ),
    m_ActiveImagesFrom( 0 ),
    m_ActiveImagesTo( 0 ),
    m_ActiveXformsFrom( 0 ),
    m_ActiveXformsTo( 0 ),
    m_TemplateNumberOfPixels( 0 ), 
    m_TemplateNumberOfSamples( 0 ), 
    m_UseTemplateData( false ),
    m_ProbabilisticSampleDensity( -1.0 ),
    m_ProbabilisticSampleUpdatesAfter( 1000000 ),
    m_ProbabilisticSampleUpdatesSince( 0 ),
    m_GaussianSmoothImagesSigma( -1.0 ),
    m_UserBackgroundValue( 0 ),
    m_UserBackgroundFlag( false ),
    m_ParametersPerXform( 0 ),
    m_RepeatIntensityHistogramMatching( false )
{
  this->m_NumberOfThreads = ThreadPool::GetGlobalThreadPool().GetNumberOfThreads();
  this->m_NumberOfTasks = 4 * this->m_NumberOfThreads - 3;

  this->m_Data.clear();
}

GroupwiseRegistrationFunctionalBase::~GroupwiseRegistrationFunctionalBase()
{
  if ( this->m_Data.size() )
    {
    const size_t numberOfImages = this->m_ImageVector.size();
    for ( size_t i = 0; i < numberOfImages; ++i )
      {
      if ( this->m_Data[i] )
	Memory::ArrayC::Delete( this->m_Data[i] );
      }
    }
}

void
GroupwiseRegistrationFunctionalBase::CreateTemplateGridFromTargets
( const std::vector<UniformVolume::SmartPtr>& targets, const int downsample )
{
  Types::Coordinate templateSize[3] = {0,0,0};
  UniformVolume::IndexType templateDims;
  Types::Coordinate templateDelta = 1e10;

  for ( size_t i = 0; i < targets.size(); ++i )
    {
    for ( int dim = 0; dim < 3; ++dim )
      {
      templateSize[dim] = std::max( templateSize[dim], targets[i]->m_Size[dim] );
      }
    templateDelta = std::min( templateDelta, targets[i]->GetMinDelta() );
    }
  
  for ( int dim = 0; dim < 3; ++dim )
    {
    templateDims[dim] = 1 + static_cast<int>( templateSize[dim] / templateDelta );
    templateSize[dim] = (templateDims[dim]-1) * templateDelta;
    }
  
  UniformVolume::SmartPtr templateGrid( new UniformVolume( templateDims, FixedVector<3,Types::Coordinate>::FromPointer( templateSize ) ) );
  this->SetTemplateGrid( templateGrid, downsample );
}

void
GroupwiseRegistrationFunctionalBase::CreateTemplateGrid
( const DataGrid::IndexType& dims, const UniformVolume::CoordinateVectorType& deltas )
{
  UniformVolume::SmartPtr templateGrid( new UniformVolume( dims, deltas ) );
  this->SetTemplateGrid( templateGrid );
}

void
GroupwiseRegistrationFunctionalBase::SetTemplateGrid
( UniformVolume::SmartPtr& templateGrid,
  const int downsample,
  const bool useTemplateData )
{ 
  this->m_TemplateGrid = templateGrid->Clone();
  this->m_UseTemplateData = useTemplateData;
  
  if ( this->m_UseTemplateData && ! this->m_TemplateGrid->GetData() )
    {
    UniformVolume::SmartPtr readImage( VolumeIO::ReadOriented( templateGrid->GetMetaInfo(META_FS_PATH) ) );
    this->m_TemplateGrid->SetData( readImage->GetData() );
    }
  
  if ( ! this->m_TemplateGrid->MetaKeyExists( META_IMAGE_ORIENTATION ) )
    {
    this->m_TemplateGrid->SetMetaInfo( META_IMAGE_ORIENTATION, AnatomicalOrientation::ORIENTATION_STANDARD );
    }
  if ( ! this->m_TemplateGrid->MetaKeyExists( META_IMAGE_ORIENTATION_ORIGINAL ) )
    {
    this->m_TemplateGrid->SetMetaInfo( META_IMAGE_ORIENTATION_ORIGINAL, AnatomicalOrientation::ORIENTATION_STANDARD );
    }
  if ( ! this->m_TemplateGrid->MetaKeyExists( META_SPACE ) )
    {
    this->m_TemplateGrid->SetMetaInfo( META_SPACE, AnatomicalOrientation::ORIENTATION_STANDARD );
    }
  if ( ! this->m_TemplateGrid->MetaKeyExists( META_SPACE_ORIGINAL ) )
    {
    this->m_TemplateGrid->SetMetaInfo( META_SPACE_ORIGINAL, AnatomicalOrientation::ORIENTATION_STANDARD );
    }

  if ( this->m_UseTemplateData )
    {
    this->m_TemplateGrid = UniformVolume::SmartPtr( this->PrepareSingleImage( this->m_TemplateGrid ) );
    }
  
  if ( downsample > 1 )
    {
    this->m_TemplateGrid = UniformVolume::SmartPtr( this->m_TemplateGrid->GetDownsampledAndAveraged( downsample, true /*approxIsotropic*/ ) );
    }
  this->m_TemplateNumberOfPixels = this->m_TemplateGrid->GetNumberOfPixels();  
  
  if ( this->m_UseTemplateData )
    {
    this->CopyTemplateData();
    }
  
  this->PrepareTargetImages();
}


void
GroupwiseRegistrationFunctionalBase
::AllocateStorage()
{
  if ( !this->m_TemplateGrid )
    {
    StdErr << "FATAL: must set template grid for groupwise registration before allocating storage\n";
    exit( 1 );
    }

  const size_t numberOfImages = this->m_OriginalImageVector.size();
  
  if ( this->m_TemplateNumberOfPixels )
    {
    if ( (this->m_ProbabilisticSampleDensity > 0) && (this->m_ProbabilisticSampleDensity < 1) )
      this->m_TemplateNumberOfSamples = static_cast<size_t>( this->m_ProbabilisticSampleDensity * this->m_TemplateNumberOfPixels);
    else
      this->m_TemplateNumberOfSamples = this->m_TemplateNumberOfPixels;
    
    if ( this->m_Data.size() )
      {
      for ( size_t i = 0; i < numberOfImages; ++i )
	{
	if ( this->m_Data[i] )
	  Memory::ArrayC::Delete( this->m_Data[i] );
	}
      }
    
    this->m_Data.resize( numberOfImages );
    for ( size_t i = 0; i < numberOfImages; ++i )
      {
      this->m_Data[i] = Memory::ArrayC::Allocate<byte>( this->m_TemplateNumberOfSamples );
      }
    
    this->m_TempData.resize( this->m_TemplateNumberOfSamples );
    }
}

void
GroupwiseRegistrationFunctionalBase
::SetTargetImages
( std::vector<UniformVolume::SmartPtr>& tImages )
{
  this->m_OriginalImageVector = tImages;

  this->m_ActiveImagesFrom = 0;
  this->m_ActiveImagesTo = tImages.size();

  this->m_ActiveXformsFrom = 0;
  this->m_ActiveXformsTo = tImages.size();

  this->m_ProbabilisticSampleUpdatesSince = 0;
}

UniformVolume::SmartPtr
GroupwiseRegistrationFunctionalBase
::PrepareSingleImage( UniformVolume::SmartPtr& image )
{
  if ( !image->GetData() )
    {
    UniformVolume::SmartPtr readImage( VolumeIO::ReadOriented( image->GetMetaInfo( META_FS_PATH ) ) );
    image->SetData( readImage->GetData() );
    }
  
  TypedArray::SmartPtr data;
  if ( this->m_GaussianSmoothImagesSigma > 0 )
    {
    data = UniformVolumeGaussianFilter( image ).GetFiltered3D( Units::GaussianSigma( this->m_GaussianSmoothImagesSigma * this->m_TemplateGrid->GetMinDelta() ) );
    
    if ( this->m_FreeAndRereadImages )
      {
      image->SetData( TypedArray::SmartPtr::Null() );
      }
    }
  else
    {
    if ( this->m_FreeAndRereadImages )
      {
      data = image->GetData();
      image->SetData( TypedArray::SmartPtr::Null() );
      }
    else
      {
      data = image->GetData()->Clone();
      }
    }
  
  UniformVolume::SmartPtr newTargetImage = image->CloneGrid();
  newTargetImage->SetData( data );
  return newTargetImage;
}

void
GroupwiseRegistrationFunctionalBase
::PrepareTargetImages()
{
  this->m_ImageVector.resize( this->m_OriginalImageVector.size() );
  for ( size_t i = 0; i < this->m_OriginalImageVector.size(); ++i )
    {
    this->m_ImageVector[i] = this->PrepareSingleImage( this->m_OriginalImageVector[i] );
    }
}

void
GroupwiseRegistrationFunctionalBase::GetParamVector( CoordinateVector& v )
{
  v.SetDim( this->ParamVectorDim() );

  for ( size_t idx = 0; idx < this->m_XformVector.size(); ++idx )
    {
    this->m_XformVector[idx]->GetParamVector( v, idx * this->m_ParametersPerXform );
    }
}

void
GroupwiseRegistrationFunctionalBase::SetParamVector( CoordinateVector& v )
{
  size_t offset = 0;
  for ( size_t xIdx = 0; xIdx < this->m_XformVector.size(); ++xIdx )
    {
    CoordinateVector vv( this->m_ParametersPerXform, v.Elements + offset, false /*free*/ );
    offset += this->m_ParametersPerXform;
    this->m_XformVector[xIdx]->SetParamVector( vv );
    }
}

void
GroupwiseRegistrationFunctionalBase::SetParamVector
( CoordinateVector& v, const size_t xformIdx )
{
  const size_t offset = this->m_ParametersPerXform * xformIdx;
  CoordinateVector vv( this->m_ParametersPerXform, v.Elements + offset, false /*free*/ );
  this->m_XformVector[xformIdx]->SetParamVector( vv );
}

void
GroupwiseRegistrationFunctionalBase
::SetParameter( const size_t param, const Types::Coordinate value )
{
  this->m_XformVector[param / this->m_ParametersPerXform]->SetParameter( param % this->m_ParametersPerXform, value );
}

void
GroupwiseRegistrationFunctionalBase
::SetParameter( const size_t xform, const size_t param, const Types::Coordinate value )
{
  this->m_XformVector[xform]->SetParameter( param, value );
}

GroupwiseRegistrationFunctionalBase::ReturnType
GroupwiseRegistrationFunctionalBase::EvaluateAt( CoordinateVector& v )
{
  if ( (this->m_ProbabilisticSampleDensity > 0) && (this->m_ProbabilisticSampleDensity < 1) )
    {
    if ( !this->m_ProbabilisticSampleUpdatesSince )
      this->UpdateProbabilisticSamples();
    (++this->m_ProbabilisticSampleUpdatesSince) %= this->m_ProbabilisticSampleUpdatesAfter;
    }

  this->SetParamVector( v );
  this->InterpolateAllImages();
  
  return this->Evaluate();
}

GroupwiseRegistrationFunctionalBase::ReturnType
GroupwiseRegistrationFunctionalBase::EvaluateWithGradient
( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step )
{
  const Self::ReturnType baseValue = this->EvaluateAt( v );

  for ( size_t param = 0; param < this->ParamVectorDim(); ++param )
    {
    g[param] = 0.0;

    const size_t imageIndex = param / this->m_ParametersPerXform;
    const size_t paramIndex = param % this->m_ParametersPerXform;

    const Types::Coordinate pStep = this->GetParamStep( param, step );
    if ( pStep > 0 )
      {
      byte* tmp = this->m_Data[imageIndex];
      this->m_Data[imageIndex] = &(this->m_TempData[0]);

      const Types::Coordinate p0 = v[param];

      this->SetParameter( imageIndex, paramIndex, p0 + pStep );
      this->InterpolateImage( imageIndex, this->m_Data[imageIndex] );
      const Self::ReturnType upper = this->Evaluate();

      this->SetParameter( imageIndex, paramIndex, p0 - pStep );
      this->InterpolateImage( imageIndex, this->m_Data[imageIndex] );
      const Self::ReturnType lower = this->Evaluate();

      this->m_Data[imageIndex] = tmp;
      this->SetParameter( imageIndex, paramIndex, p0 );

      if ( (upper > baseValue) || (lower > baseValue) )
	{
	g[param] = (upper - lower);
	}
      }
    }

  if ( this->m_ForceZeroSum )
    {
    this->ForceZeroSumGradient( g );
    }

  return baseValue;
}

void
GroupwiseRegistrationFunctionalBase::UpdateProbabilisticSamples()
{
  this->m_ProbabilisticSamples.resize( this->m_TemplateNumberOfSamples );
  
  const size_t firstSample = 0;
  const size_t lastSample = this->m_TemplateNumberOfSamples;

  for ( size_t i = firstSample; i < lastSample; ++i )
    {
    const size_t sample = static_cast<size_t>( this->m_TemplateNumberOfPixels * MathUtil::UniformRandom() );
    this->m_ProbabilisticSamples[i] = sample;
    }
}

void
GroupwiseRegistrationFunctionalBase
::InterpolateAllImages()
{
  for ( size_t idx = this->m_ActiveImagesFrom; idx < this->m_ActiveImagesTo; ++idx )
    {
    this->InterpolateImage( idx, this->m_Data[idx] );
    }
}

void
GroupwiseRegistrationFunctionalBase
::ForceZeroSumGradient( CoordinateVector& g ) const
{
  const size_t numberOfXforms = this->m_XformVector.size();
  const size_t zeroSumFirstN = this->m_ForceZeroSumFirstN ? this->m_ForceZeroSumFirstN : numberOfXforms;

#pragma omp parallel for
  for ( int param = 0; param < static_cast<int>( this->m_ParametersPerXform ); ++param )
    {
    Types::Coordinate avg = 0;

    size_t pparam = param;
    for ( size_t idx = 0; idx < zeroSumFirstN; ++idx, pparam += this->m_ParametersPerXform )
      {
      avg += g[pparam];
      }
    
    avg *= 1.0 / zeroSumFirstN;
    
    pparam = param;
    for ( size_t idx = 0; idx < numberOfXforms; ++idx, pparam += this->m_ParametersPerXform  )
      {
      g[pparam] -= avg;
      }
    }

  // if the remaining vector is smaller than threshold, throw it out.
  const Types::Coordinate gMax = g.MaxNorm();
  if ( gMax < 1e-3 )
    {
    g.Clear();
    }
}

bool
GroupwiseRegistrationFunctionalBase
::Wiggle()
{
  bool wiggle = false;

  if ( (this->m_ProbabilisticSampleDensity > 0) && (this->m_ProbabilisticSampleDensity < 1) )
    {
    this->m_ProbabilisticSampleUpdatesSince = 0;
    wiggle = true;
    }

  if ( this->m_RepeatIntensityHistogramMatching )
    {
    TypedArray::SmartPtr referenceData = this->m_TemplateGrid->GetData();
    if ( !this->m_UseTemplateData )
      referenceData = TypedArray::SmartPtr::Null();
    
    for ( size_t i = 0; i < this->m_OriginalImageVector.size(); ++i )
      {
      UniformVolume::SmartPtr scaledImage;
      if ( this->m_OriginalImageVector[i]->GetData() )
	{
	scaledImage = UniformVolume::SmartPtr( this->m_OriginalImageVector[i]->Clone( true /*copyData*/ ) );
	}
      else
	{
	scaledImage = UniformVolume::SmartPtr( VolumeIO::ReadOriented( this->m_OriginalImageVector[i]->GetMetaInfo( META_FS_PATH ) ) );
	}

      UniformVolume::SmartPtr reformatImage( this->GetReformattedImage( scaledImage, i ) );
      if ( referenceData )
	{
	scaledImage->GetData()->ApplyFunctionObject( TypedArrayFunctionHistogramMatching( *(reformatImage->GetData()), *referenceData ) );
	}
      else
	{
	referenceData = reformatImage->GetData();
	}
      
      this->m_ImageVector[i] = UniformVolume::SmartPtr( this->PrepareSingleImage( scaledImage ) );
      }
    this->InterpolateAllImages();
    wiggle = true;
    }

  return wiggle;
}

void
GroupwiseRegistrationFunctionalBase
::CopyTemplateData()
{
  const TypedArray* dataArray = this->m_TemplateGrid->GetData();

  if ( dataArray )
    {
    const size_t size = dataArray->GetDataSize();
    this->m_TemplateData.resize( size );
    
    for ( size_t i = 0; i < size; ++i )
      {
      Types::DataItem value;
      if ( dataArray->Get( value, i ) )
	this->m_TemplateData[i] = static_cast<byte>( value );
      else
	this->m_TemplateData[i] = this->m_PaddingValue;
      }
    }
}

void
GroupwiseRegistrationFunctionalBase
::DebugWriteImages()
{
  this->InterpolateAllImages();
  UniformVolume::SmartPtr writeVolume( this->m_TemplateGrid->CloneGrid() );
  writeVolume->CreateDataArray( TYPE_BYTE );

  for ( size_t i = 0; i < this->m_TemplateNumberOfPixels; ++i )
    {
    writeVolume->SetDataAt( this->m_TemplateData[i], i );
    }
  VolumeIO::Write( *writeVolume, "template.nii" );

  for ( size_t n = 0; n < this->m_ImageVector.size(); ++n )
    {
    for ( size_t i = 0; i < this->m_TemplateNumberOfPixels; ++i )
      {
      writeVolume->SetDataAt( this->m_Data[n][i], i );
      }

    char path[PATH_MAX];
    sprintf( path, "target%02d.nii", static_cast<int>( n ) );
    VolumeIO::Write( *writeVolume, path );
    }
}

const UniformVolume::SmartPtr
GroupwiseRegistrationFunctionalBase
::GetReformattedImage( const UniformVolume::SmartPtr& targetGrid, const size_t idx ) const
{
  ReformatVolume reformat;
  reformat.SetInterpolation( Interpolators::LINEAR );
  reformat.SetReferenceVolume( targetGrid );
  reformat.SetFloatingVolume( this->m_OriginalImageVector[idx] );
  
  reformat.SetWarpXform( WarpXform::SmartPtr::DynamicCastFrom( this->m_XformVector[idx] ) );
  reformat.SetAffineXform( AffineXform::SmartPtr::DynamicCastFrom( this->m_XformVector[idx] ) );
  
  if ( this->m_UserBackgroundFlag )
    {
    reformat.SetPaddingValue( this->m_UserBackgroundValue );
    }
  
  UniformVolume::SmartPtr result = reformat.PlainReformat();

  if ( this->m_UserBackgroundFlag )
    {
    result->GetData()->ClearPaddingFlag();
    }
  return result;
}

} // namespace cmtk
