/*
//
//  Copyright 2016 Google, Inc.
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
//  $Revision: 1768 $
//
//  $LastChangedDate: 2010-05-27 15:16:01 -0700 (Thu, 27 May 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <Registration/cmtkGroupwiseRegistrationFunctionalXformTemplate.h>

#include <System/cmtkDebugOutput.h>

#include <Base/cmtkRegionIndexIterator.h>

namespace
cmtk
{

GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::GroupwiseRegistrationFunctionalXformTemplate()
  : m_MaximumNumberOfPixelsVOI( 0 ), m_MaximumNumberOfPixelsPerLineVOI( 0 ), m_ForceZeroSumNoAffine( false )
{
  this->m_ParametersPerXform = 0;
  this->m_WarpFastMode = true;
  this->m_PartialGradientMode = false;
  this->m_PartialGradientThreshold = static_cast<Types::DataItem>( 0.01 );
  this->m_DeactivateUninformativeMode = false;
  this->m_NumberOfActiveControlPoints = 0;
  this->m_JacobianConstraintWeight = 0.0;
  this->m_BendingEnergyWeight = 0.0;

  this->SetFreeAndRereadImages( true );
}

void
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::SetTemplateGrid
( UniformVolume::SmartPtr& templateGrid, const int downsample, const bool useTemplateData )
{
  this->Superclass::SetTemplateGrid( templateGrid, downsample, useTemplateData );
  
  if ( this->m_XformVector.size() )
    {
    for ( size_t i = 0; i < this->m_XformVector.size(); ++i )
      {
      dynamic_cast<SplineWarpXform&>( *(this->m_XformVector[i]) ).RegisterVolume( *(this->m_TemplateGrid) );	
      }      
    this->UpdateVolumesOfInfluence();
    }
}

void
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::InitializeXformsFromAffine
( const Types::Coordinate gridSpacing, std::vector<AffineXform::SmartPtr> initialAffineXformsVector, const bool exactSpacing )
{
  this->m_InitialAffineXformsVector = initialAffineXformsVector;  

  this->m_XformVector.resize( this->m_ImageVector.size() );
  this->m_InitialRotationsVector.resize( this->m_ImageVector.size() );

  for ( size_t i = 0; i < this->m_ImageVector.size(); ++i )
    {
    SplineWarpXform::SmartPtr xform( new SplineWarpXform( this->m_TemplateGrid->m_Size, gridSpacing, initialAffineXformsVector[i], exactSpacing ) );
    xform->RegisterVolume( *(this->m_TemplateGrid) );
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
}

void
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::RefineTransformationGrids()
{
  for ( size_t i = 0; i < this->m_XformVector.size(); ++i )
    {
    this->GetXformByIndex(i)->Refine();
    dynamic_cast<SplineWarpXform&>( *(this->m_XformVector[i]) ).RegisterVolume( *(this->m_TemplateGrid) );
    }

  this->m_ParametersPerXform = this->m_XformVector[0]->VariableParamVectorDim();

  this->UpdateVolumesOfInfluence();
}

void
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::UpdateVolumesOfInfluence()
{
  const UniformVolume::CoordinateRegionType templateDomain( this->m_TemplateGrid->m_Offset, this->m_TemplateGrid->m_Offset + this->m_TemplateGrid->m_Size );
  
  this->m_VolumeOfInfluenceArray.resize( this->m_ParametersPerXform / 3 );
  
  this->m_MaximumNumberOfPixelsPerLineVOI = 0;
  this->m_MaximumNumberOfPixelsVOI = 0;
  
  const SplineWarpXform& xform0 = *(this->GetXformByIndex(0));
  for ( size_t param = 0; param < this->m_ParametersPerXform; param += 3 ) 
    { 
    DataGrid::RegionType& voi = this->m_VolumeOfInfluenceArray[param/3];
    voi = this->m_TemplateGrid->GetGridRange( xform0.GetVolumeOfInfluence( param, templateDomain ) );
    
    this->m_MaximumNumberOfPixelsVOI = std::max<size_t>( voi.Size(), this->m_MaximumNumberOfPixelsVOI );
    this->m_MaximumNumberOfPixelsPerLineVOI = std::max<size_t>( voi.To()[0]-voi.From()[0], this->m_MaximumNumberOfPixelsPerLineVOI );
    }
}
  
bool
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::UpdateParamStepArray()
{
  bool changed = false;

  this->m_ParamStepArray.resize( this->ParamVectorDim() );

  if ( ( this->m_DeactivateUninformativeMode || this->m_DisableControlPointsMask ) &&
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
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
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
	for ( int param = 0; param < static_cast<int>( this->m_ParametersPerXform ); param += 3 )
	  {
	  const FixedVector<3,Types::Coordinate> projected( affineXform->RotateScaleShear( FixedVector<3,Types::Coordinate>::FromPointer( gX+param ) ) );
	  for ( size_t i = 0; i<3; ++i )
	    gX[param+i] = projected[i];
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
	for ( int param = 0; param < static_cast<int>( this->m_ParametersPerXform ); param += 3 )
	  {
	  const FixedVector<3,Types::Coordinate> projected( affineXform->RotateScaleShear( FixedVector<3,Types::Coordinate>::FromPointer( gX+param ) ) );
	  for ( size_t i = 0; i<3; ++i )
	    gX[param+i] = projected[i];
	  }
	}
      }
    }
}

bool
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::Wiggle()
{
  bool wiggle = this->Superclass::Wiggle();

  if ( this->m_PartialGradientMode )
    {
    wiggle = wiggle || this->UpdateParamStepArray();
    }

  return wiggle;
}

void
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::InterpolateImage
( const size_t idx, byte* const destination )
{
  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  const size_t numberOfThreads = threadPool.GetNumberOfThreads();
  std::vector<InterpolateImageThreadParameters> params( numberOfThreads );

  for ( size_t thread = 0; thread < numberOfThreads; ++thread )
    {
    params[thread].thisObject = this;    
    params[thread].m_Idx = idx;    
    params[thread].m_Destination = destination;    
    }
  
  threadPool.Run( InterpolateImageThread, params );
}

void
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::InterpolateImageThread
( void* args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t )
{
  InterpolateImageThreadParameters* threadParameters = static_cast<InterpolateImageThreadParameters*>( args );
  
  const Self* This = threadParameters->thisObject;
  const size_t idx = threadParameters->m_Idx;
  byte* destination = threadParameters->m_Destination;

  const SplineWarpXform* xform = This->GetXformByIndex(idx);
  const UniformVolume* target = This->m_ImageVector[idx];
  const byte* dataPtr = static_cast<const byte*>( target->GetData()->GetDataPtr() );

  const byte paddingValue = This->m_PaddingValue;
  const byte backgroundValue = This->m_UserBackgroundFlag ? This->m_PrivateUserBackgroundValue : paddingValue;

  const Types::GridIndexType dimsX = This->m_TemplateGrid->GetDims()[AXIS_X];
  const Types::GridIndexType dimsY = This->m_TemplateGrid->GetDims()[AXIS_Y];
  const Types::GridIndexType dimsZ = This->m_TemplateGrid->GetDims()[AXIS_Z];

  std::vector<Xform::SpaceVectorType> v( dimsX );
  byte value;

  const Types::GridIndexType rowCount = ( dimsY * dimsZ );
  const Types::GridIndexType rowFrom = ( rowCount / taskCnt ) * taskIdx;
  const Types::GridIndexType rowTo = ( taskIdx == (taskCnt-1) ) ? rowCount : ( rowCount / taskCnt ) * ( taskIdx + 1 );
  Types::GridIndexType rowsToDo = rowTo - rowFrom;
  
  Types::GridIndexType yFrom = rowFrom % dimsY;
  Types::GridIndexType zFrom = rowFrom / dimsY;
  
  byte *wptr = destination + rowFrom * dimsX;
  for ( Types::GridIndexType z = zFrom; (z < dimsZ) && rowsToDo; ++z ) 
    {
    for ( Types::GridIndexType y = yFrom; (y < dimsY) && rowsToDo; yFrom = 0, ++y, --rowsToDo )
      {
      xform->GetTransformedGridRow( dimsX, &v[0], 0, y, z );
      for ( Types::GridIndexType x = 0; x < dimsX; ++x )
	{
	if ( target->ProbeData( value, dataPtr, v[x] ) )
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
}

void
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>::UpdateActiveControlPoints()
{
  const size_t numberOfControlPoints = this->m_VolumeOfInfluenceArray.size();
  
  if ( numberOfControlPoints )
    {
    this->m_ActiveControlPointFlags.resize( numberOfControlPoints );
    std::fill( this->m_ActiveControlPointFlags.begin(), this->m_ActiveControlPointFlags.end(), true );
    this->m_NumberOfActiveControlPoints = numberOfControlPoints;
    }

  if ( this->m_DisableControlPointsMask )
    {
    size_t cntDisabled = 0;

    const UniformVolume::CoordinateRegionType templateDomain( this->m_TemplateGrid->m_Offset, this->m_TemplateGrid->m_Offset + this->m_TemplateGrid->m_Size );
    const SplineWarpXform& xform0 = *(this->GetXformByIndex(0));
    for ( size_t cp = 0; cp < numberOfControlPoints; ++cp )
      {
      const DataGrid::RegionType maskRegion = this->m_DisableControlPointsMask->GetGridRange( xform0.GetVolumeOfInfluence( 3*cp, templateDomain, false /*force slow, accurate mode*/ ) );
      for ( RegionIndexIterator<DataGrid::RegionType> it( maskRegion ); it != it.end(); ++it )
	{
	if ( this->m_DisableControlPointsMask->GetDataAt( this->m_DisableControlPointsMask->GetOffsetFromIndex( it.Index() ) ) > 0 )
	  {
	  this->m_ActiveControlPointFlags[cp] = false;
	  ++cntDisabled;
	  break;
	  }
	}
      }

    DebugOutput( 2 ) << "Disabled " << cntDisabled << " control points due to provided mask.\n";
    }
}


} // namespace cmtk
