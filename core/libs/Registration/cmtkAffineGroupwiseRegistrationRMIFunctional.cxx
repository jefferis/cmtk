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

#include <cmtkAffineGroupwiseRegistrationRMIFunctional.h>

#include <cmtkMathUtil.h>
#include <cmtkVolumeAxesHash.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

AffineGroupwiseRegistrationRMIFunctional::AffineGroupwiseRegistrationRMIFunctional()
  : m_XformNumberDOFs( 9 )
{
  this->m_ParametersPerXform = AffineXform::TotalNumberOfParameters;
}

AffineGroupwiseRegistrationRMIFunctional::~AffineGroupwiseRegistrationRMIFunctional()
{
}

void
AffineGroupwiseRegistrationRMIFunctional
::SetXformNumberDOFs( const int numberDOFs )
{ 
  this->m_XformNumberDOFs = numberDOFs;
  std::vector<Xform::SmartPtr>::iterator it = this->m_XformVector.begin();
  while ( it != this->m_XformVector.end() )
    {
    AffineXform::SmartPtr::DynamicCastFrom(*it)->SetNumberDOFs( this->m_XformNumberDOFs );
    ++it;
    }
}

void
AffineGroupwiseRegistrationRMIFunctional::InitializeXforms
( const bool alignCenters )
{
  const Vector3D centerTemplate = this->m_TemplateGrid->GetCenterCropRegion();
  
  this->m_XformVector.resize( this->m_ImageVector.size() );
  for ( size_t i = 0; i < this->m_ImageVector.size(); ++i )
    {
    AffineXform::SmartPtr xform( new AffineXform );
    xform->SetNumberDOFs( this->m_XformNumberDOFs );
    xform->SetUseLogScaleFactors( true );
    xform->SetCenter( centerTemplate.XYZ );
    this->m_XformVector[i] = xform;
 
    if ( alignCenters )
      {
      Vector3D center = this->m_ImageVector[i]->GetCenter();
      center -= centerTemplate;
      this->GetXformByIndex(i)->SetXlate( center.XYZ );
      }
    }
}

void
AffineGroupwiseRegistrationRMIFunctional::SetXforms
( const std::vector<AffineXform::SmartPtr>& xformVector )
{
  const Vector3D centerTemplate = this->m_TemplateGrid->GetCenterCropRegion();
  
  this->m_XformVector.resize( xformVector.size() );
  for ( size_t i = 0; i < this->m_XformVector.size(); ++i )
    {
    AffineXform::SmartPtr xform( new AffineXform( *(xformVector[i]) ) );
    xform->SetNumberDOFs( this->m_XformNumberDOFs );
    xform->SetUseLogScaleFactors( true );

//    const Vector3D center = this->m_ImageVector[i]->GetCenterCropRegion();
    xform->ChangeCenter( centerTemplate.XYZ );

    this->m_XformVector[i] = xform;
    }
}

void
AffineGroupwiseRegistrationRMIFunctional::InterpolateImage
( const size_t idx, byte* const destination )
{
  const VolumeAxesHash gridHash( *this->m_TemplateGrid, this->GetXformByIndex(idx)->GetInverse() );

  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  ThreadParameterArray<Self,InterpolateImageThreadParameters> params( this, numberOfThreads );

  for ( size_t thread = 0; thread < numberOfThreads; ++thread )
    {
    params[thread].m_Idx = idx;    
    params[thread].m_Destination = destination;    
    params[thread].m_HashX = gridHash[0];
    params[thread].m_HashY = gridHash[1];
    params[thread].m_HashZ = gridHash[2];
    }
  
  if ( (this->m_ProbabilisticSampleDensity > 0) && (this->m_ProbabilisticSampleDensity < 1) )
    params.RunInParallel( &InterpolateImageProbabilisticThread );
  else
    params.RunInParallel( &InterpolateImageThread );    
}

CMTK_THREAD_RETURN_TYPE
AffineGroupwiseRegistrationRMIFunctional::InterpolateImageThread
( void* args )
{
  InterpolateImageThreadParameters* threadParameters = static_cast<InterpolateImageThreadParameters*>( args );
  
  const Self* This = threadParameters->thisObject;
  const int threadID = threadParameters->ThisThreadIndex;
  const int numberOfThreads = threadParameters->NumberOfThreads;
  const size_t idx = threadParameters->m_Idx;
  byte* destination = threadParameters->m_Destination;

  const UniformVolume* target = This->m_ImageVector[idx];

  const byte paddingValue = This->m_PaddingValue;
  const byte backgroundValue = This->m_UserBackgroundFlag ? This->m_UserBackgroundValue : paddingValue;

  Vector3D v;
  byte value;
  const byte* dataPtr = static_cast<const byte*>( target->GetData()->GetDataPtr() );

  const int dimsX = This->m_TemplateGrid->GetDims( AXIS_X );
  const int dimsY = This->m_TemplateGrid->GetDims( AXIS_Y );
  const int dimsZ = This->m_TemplateGrid->GetDims( AXIS_Z );

  const int rowCount = ( dimsY * dimsZ );
  const int rowFrom = ( rowCount / numberOfThreads ) * threadID;
  const int rowTo = ( threadID == (numberOfThreads-1) ) ? rowCount : ( rowCount / numberOfThreads ) * ( threadID + 1 );
  int rowsToDo = rowTo - rowFrom;
  
  int yFrom = rowFrom % dimsY;
  int zFrom = rowFrom / dimsY;
  
  Vector3D planeStart, rowStart;
  byte *wptr = destination + rowFrom * dimsX;
  for ( int z = zFrom; (z < dimsZ) && rowsToDo; ++z ) 
    {
    planeStart = threadParameters->m_HashZ[z];
    for ( int y = yFrom; (y < dimsY) && rowsToDo; yFrom = 0, ++y, --rowsToDo )
      {
      (rowStart = planeStart) += threadParameters->m_HashY[y];
      for ( int x = 0; x < dimsX; ++x )
	{
	(v = rowStart) += threadParameters->m_HashX[x];	
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

CMTK_THREAD_RETURN_TYPE
AffineGroupwiseRegistrationRMIFunctional::InterpolateImageProbabilisticThread
( void* args )
{
  InterpolateImageThreadParameters* threadParameters = static_cast<InterpolateImageThreadParameters*>( args );
  
  const Self* This = threadParameters->thisObject;
  const int threadID = threadParameters->ThisThreadIndex;
  const int numberOfThreads = threadParameters->NumberOfThreads;
  const size_t idx = threadParameters->m_Idx;
  byte* destination = threadParameters->m_Destination;

  const AffineXform* xform = This->GetXformByIndex(idx);
  const UniformVolume* target = This->m_ImageVector[idx];

  const byte paddingValue = This->m_PaddingValue;
  const byte backgroundValue = This->m_UserBackgroundFlag ? This->m_UserBackgroundValue : paddingValue;

  Vector3D v;
  byte value;
  const byte* dataPtr = static_cast<const byte*>( target->GetData()->GetDataPtr() );

  const size_t samplesPerThread = This->m_ProbabilisticSamples.size() / numberOfThreads;
  const size_t startIdx = threadID * samplesPerThread;
  const size_t endIdx = std::min( startIdx + samplesPerThread, This->m_ProbabilisticSamples.size() );

  byte *wptr = destination + startIdx;
  for ( size_t i = startIdx; i < endIdx; ++i, ++wptr )
    {
    const size_t offset = This->m_ProbabilisticSamples[i];
    This->m_TemplateGrid->GetGridLocation( v, offset );
    xform->ApplyInPlace( v );
    
    if ( target->ProbeData( value, dataPtr, v ) )
      {
      *wptr = value;
      }
    else
      {
      *wptr = backgroundValue;
      }
    }

  return CMTK_THREAD_RETURN_VALUE;
}

} // namespace cmtk
