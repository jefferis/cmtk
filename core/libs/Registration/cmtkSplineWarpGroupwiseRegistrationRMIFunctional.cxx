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

#include <cmtkSplineWarpGroupwiseRegistrationRMIFunctional.h>

#include <cmtkMathUtil.h>
#include <cmtkMatrix.h>

#include <cmtkThreadParameterArray.h>

#include <algorithm>

#ifdef CMTK_BUILD_MPI
#  include <mpi.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

SplineWarpGroupwiseRegistrationRMIFunctional
::SplineWarpGroupwiseRegistrationRMIFunctional()
{
  this->m_ParametersPerXform = 0;
  this->m_WarpFastMode = true;
  this->m_PartialGradientMode = false;
  this->m_PartialGradientThreshold = 0.01;
  this->m_DeactivateUninformativeMode = false;
  this->m_NumberOfActiveControlPoints = 0;
  this->m_NeedsUpdateInformationByControlPoint = true;

  this->SetFreeAndRereadImages( true );
}

SplineWarpGroupwiseRegistrationRMIFunctional
::~SplineWarpGroupwiseRegistrationRMIFunctional()
{
}

void
SplineWarpGroupwiseRegistrationRMIFunctional
::SetTemplateGrid
( UniformVolume::SmartPtr& templateGrid, 
  const int downsample,
  const bool useTemplateData )
{
  this->Superclass::SetTemplateGrid( templateGrid, downsample, useTemplateData );
  
  if ( this->m_XformVector.size() )
    {
    for ( size_t i = 0; i < this->m_XformVector.size(); ++i )
      {
      dynamic_cast<SplineWarpXform&>( *(this->m_XformVector[i]) ).RegisterVolume( this->m_TemplateGrid );	
      }      
    this->UpdateVolumesOfInfluence();
    }
}

void
SplineWarpGroupwiseRegistrationRMIFunctional::InitializeXforms
( const Types::Coordinate gridSpacing, std::vector<AffineXform::SmartPtr> initialAffineXformsVector )
{
  this->m_InitialAffineXformsVector = initialAffineXformsVector;

  this->m_XformVector.resize( this->m_ImageVector.size() );
  for ( size_t i = 0; i < this->m_ImageVector.size(); ++i )
    {
    SplineWarpXform::SmartPtr xform
      ( new SplineWarpXform( this->m_TemplateGrid->Size, gridSpacing,
				initialAffineXformsVector[i], true /*exactSpacing*/ ) );
    xform->RegisterVolume( this->m_TemplateGrid );
    this->m_XformVector[i] = xform;
    }

  this->m_ParametersPerXform = this->m_XformVector[0]->VariableParamVectorDim();

  this->UpdateVolumesOfInfluence();
}

void
SplineWarpGroupwiseRegistrationRMIFunctional::RefineTransformationGrids()
{

  for ( size_t i = 0; i < this->m_XformVector.size(); ++i )
    {
    this->GetXformByIndex(i)->Refine();
    dynamic_cast<SplineWarpXform&>( *(this->m_XformVector[i]) ).RegisterVolume( this->m_TemplateGrid );
    }

  this->m_ParametersPerXform = this->m_XformVector[0]->VariableParamVectorDim();

  this->UpdateVolumesOfInfluence();
}

void
SplineWarpGroupwiseRegistrationRMIFunctional::UpdateProbabilisticSamples()
{
  this->Superclass::UpdateProbabilisticSamples();
}

bool
SplineWarpGroupwiseRegistrationRMIFunctional::UpdateParamStepArray()
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
SplineWarpGroupwiseRegistrationRMIFunctional
::UpdateVolumesOfInfluence()
{
  const Vector3D templateFrom( this->m_TemplateGrid->m_Offset );
  const Vector3D templateTo(  this->m_TemplateGrid->m_Offset + this->m_TemplateGrid->Size );
  
  this->m_VolumeOfInfluenceArray.resize( this->m_ParametersPerXform / 3 );

  const SplineWarpXform* xform0 = this->GetXformByIndex(0);
  for ( size_t param = 0; param < this->m_ParametersPerXform; param += 3 ) 
    { 
    Vector3D fromVOI, toVOI;
    xform0->GetVolumeOfInfluence( param, templateFrom, templateTo, fromVOI, toVOI );
    this->m_VolumeOfInfluenceArray[param/3] = this->m_TemplateGrid->GetGridRange( fromVOI, toVOI );
    }
  
  this->m_NeedsUpdateInformationByControlPoint = true;
}

void
SplineWarpGroupwiseRegistrationRMIFunctional
::UpdateInformationByControlPoint()
{
  this->m_NeedsUpdateInformationByControlPoint = false;

  const size_t numberOfControlPoints = this->m_VolumeOfInfluenceArray.size();
#ifdef CMTK_BUILD_MPI
  const size_t cpPerNode = 1 + numberOfControlPoints / this->m_SizeMPI;
  std::vector<byte> tmpInformationByCP( cpPerNode );
#else
  this->m_InformationByControlPoint.resize( numberOfControlPoints );
#endif

  const byte paddingValue = this->m_PaddingValue;

#ifdef CMTK_BUILD_MPI
  const size_t beginCP = this->m_RankMPI * cpPerNode;
  const size_t endCP = std::min<size_t>( beginCP + cpPerNode, numberOfControlPoints );
#else
  const size_t beginCP = 0;
  const size_t endCP = numberOfControlPoints;
#endif
  for ( size_t cp = beginCP; cp < endCP; ++cp ) 
    {
#ifdef CMTK_BUILD_MPI
    const size_t ofs = cp-beginCP;
    tmpInformationByCP[ofs] = 0;
#else
    this->m_InformationByControlPoint[cp] = 0;
#endif

    std::vector<DataGrid::RegionType>::const_iterator voi = this->m_VolumeOfInfluenceArray.begin() + cp;
    for ( size_t img = this->m_ActiveImagesFrom; img < this->m_ActiveImagesTo; ++img )
      {
      const byte* dataPtrImg = this->m_Data[img];
      
      byte voiMin = 255, voiMax = 0;
      for ( int z = voi->From()[2]; z < voi->To()[2]; ++z ) 
	{
	for ( int y = voi->From()[1]; y < voi->To()[1]; ++y )
	  {
	  size_t ofs = this->m_TemplateGrid->GetOffsetFromIndex( voi->From()[0], y, z );
	  for ( int x = voi->From()[0]; x < voi->To()[0]; ++x, ++ofs )
	    {
	    const byte data = dataPtrImg[ofs];
	    if ( data != paddingValue )
	      {
	      voiMin = std::min( data, voiMin );
	      voiMax = std::max( data, voiMax );
	      }
	    }
	  }
	}
#ifdef CMTK_BUILD_MPI
      tmpInformationByCP[ofs] = std::max( (byte)(voiMax-voiMin), tmpInformationByCP[ofs] );
#else
      this->m_InformationByControlPoint[cp] = std::max( (byte)(voiMax-voiMin), this->m_InformationByControlPoint[cp] );    
#endif
      }
    }

#ifdef CMTK_BUILD_MPI
  this->m_InformationByControlPoint.resize( cpPerNode * this->m_SizeMPI );
  MPI::COMM_WORLD.Allgather( &tmpInformationByCP[0], cpPerNode, MPI::CHAR, &this->m_InformationByControlPoint[0], cpPerNode, MPI::CHAR );
#endif
  
  this->UpdateActiveControlPoints();
}
  
void
SplineWarpGroupwiseRegistrationRMIFunctional::UpdateControlPointSchedule()
{
  const SplineWarpXform* xform0 = this->GetXformByIndex(0);
  this->m_ControlPointSchedule.resize( xform0->GetNumberOfControlPoints() );
  this->m_ControlPointScheduleOverlapFreeMaxLength = (xform0->m_Dims[0] / 4) * (xform0->m_Dims[1] / 4) * (xform0->m_Dims[2] / 4);
  
  size_t ofs = 0;
  for ( int z = 0; z < 4; ++z )
    {
    for ( int y = 0; y < 4; ++y )
      {
      for ( int x = 0; x < 4; ++x )
	{
	for ( int k = z; k < xform0->m_Dims[2]; k += 4 )
	  {
	  for ( int j = y; j < xform0->m_Dims[1]; j += 4 )
	    {
	    for ( int i = x; i < xform0->m_Dims[0]; i += 4, ++ofs )
	      {
	      this->m_ControlPointSchedule[ofs] = i + xform0->m_Dims[0] * ( j + xform0->m_Dims[1] * k );
	      }
	    }
	  }
	}
      }
    }
}

void
SplineWarpGroupwiseRegistrationRMIFunctional::UpdateActiveControlPoints()
{
  if ( this->m_DeactivateUninformativeMode )
    {
    const size_t numberOfControlPoints = this->m_VolumeOfInfluenceArray.size();
    
    if ( numberOfControlPoints )
      {
      this->m_ActiveControlPointFlags.resize( numberOfControlPoints );
      this->m_NumberOfActiveControlPoints = 0;
      
      const Vector3D templateFrom( this->m_TemplateGrid->m_Offset );
      const Vector3D templateTo(  this->m_TemplateGrid->m_Offset + this->m_TemplateGrid->Size );
      Vector3D fromVOI, toVOI;
      
      std::vector<DataGrid::RegionType>::const_iterator voi = this->m_VolumeOfInfluenceArray.begin();
      for ( size_t cp = 0; cp < numberOfControlPoints; ++cp, ++voi )
	{
	this->m_ActiveControlPointFlags[cp] = (this->m_InformationByControlPoint[cp] > (this->m_HistogramBins / 4) );
	if ( this->m_ActiveControlPointFlags[cp] ) 
	  ++this->m_NumberOfActiveControlPoints;
	}
      
      StdErr << "Enabled " << this->m_NumberOfActiveControlPoints 
		<< "/" << this->m_ParametersPerXform / 3
		<< " control points.\n";
      }
    }
  else
    {
    this->m_NumberOfActiveControlPoints = this->m_VolumeOfInfluenceArray.size();
    }
  
  this->UpdateParamStepArray();
  this->UpdateControlPointSchedule();
}

SplineWarpGroupwiseRegistrationRMIFunctional::ReturnType
SplineWarpGroupwiseRegistrationRMIFunctional::Evaluate()
{
  return this->Superclass::Evaluate();
}

void
SplineWarpGroupwiseRegistrationRMIFunctional
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
SplineWarpGroupwiseRegistrationRMIFunctional::InterpolateImageThread
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

  const int dimsX = This->m_TemplateGrid->GetDims()[AXIS_X];
  const int dimsY = This->m_TemplateGrid->GetDims()[AXIS_Y];
  const int dimsZ = This->m_TemplateGrid->GetDims()[AXIS_Z];

  std::vector<Vector3D> vectorList( dimsX );
  byte value;

  const byte paddingValue = This->m_PaddingValue;
  const byte backgroundValue = This->m_UserBackgroundFlag ? This->m_UserBackgroundValue : paddingValue;

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
      xform->GetTransformedGridSequence( &(vectorList[0]), dimsX, 0, y, z );
      for ( int x = 0; x < dimsX; ++x )
	{
	if ( target->ProbeData( value, dataPtr, vectorList[x] ) )
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

bool
SplineWarpGroupwiseRegistrationRMIFunctional::Wiggle()
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
#  include "cmtkSplineWarpGroupwiseRegistrationRMIFunctionalMPI.txx"
#else
#  include "cmtkSplineWarpGroupwiseRegistrationRMIFunctional.txx"
#endif
