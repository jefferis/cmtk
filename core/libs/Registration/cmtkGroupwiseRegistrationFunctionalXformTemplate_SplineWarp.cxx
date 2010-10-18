#include <Registration/cmtkGroupwiseRegistrationFunctionalXformTemplate.h>

namespace
cmtk
{

GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::GroupwiseRegistrationFunctionalXformTemplate()
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
      dynamic_cast<SplineWarpXform&>( *(this->m_XformVector[i]) ).RegisterVolume( this->m_TemplateGrid );	
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
}

void
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::RefineTransformationGrids()
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
GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>
::UpdateVolumesOfInfluence()
{
  const Vector3D templateFrom( this->m_TemplateGrid->m_Offset );
  const Vector3D templateTo(  this->m_TemplateGrid->m_Offset + this->m_TemplateGrid->Size );
  
  this->m_VolumeOfInfluenceArray.resize( this->m_ParametersPerXform / 3 );

  this->m_MaximumNumberOfPixelsPerLineVOI = 0;
  this->m_MaximumNumberOfPixelsVOI = 0;
  
  const SplineWarpXform* xform0 = this->GetXformByIndex(0);
  for ( size_t param = 0; param < this->m_ParametersPerXform; param += 3 ) 
    { 
    Vector3D fromVOI, toVOI;
    xform0->GetVolumeOfInfluence( param, templateFrom, templateTo, fromVOI, toVOI );

    DataGrid::RegionType& voi = this->m_VolumeOfInfluenceArray[param/3];
    voi = this->m_TemplateGrid->GetGridRange( fromVOI, toVOI );
    
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
	for ( size_t param = 0; param < this->m_ParametersPerXform; param += 3 )
	  {
	  const FixedVector<3,Types::Coordinate> projected( affineXform->RotateScaleShear( FixedVector<3,Types::Coordinate>( gX+param ) ) );
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
	for ( size_t param = 0; param < this->m_ParametersPerXform; param += 3 )
	  {
	  const FixedVector<3,Types::Coordinate> projected( affineXform->RotateScaleShear( FixedVector<3,Types::Coordinate>( gX+param ) ) );
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

  const int dimsX = This->m_TemplateGrid->GetDims()[AXIS_X];
  const int dimsY = This->m_TemplateGrid->GetDims()[AXIS_Y];
  const int dimsZ = This->m_TemplateGrid->GetDims()[AXIS_Z];

  std::vector<Xform::SpaceVectorType> v( dimsX );
  byte value;

  const int rowCount = ( dimsY * dimsZ );
  const int rowFrom = ( rowCount / taskCnt ) * taskIdx;
  const int rowTo = ( taskIdx == (taskCnt-1) ) ? rowCount : ( rowCount / taskCnt ) * ( taskIdx + 1 );
  int rowsToDo = rowTo - rowFrom;
  
  int yFrom = rowFrom % dimsY;
  int zFrom = rowFrom / dimsY;
  
  byte *wptr = destination + rowFrom * dimsX;
  for ( int z = zFrom; (z < dimsZ) && rowsToDo; ++z ) 
    {
    for ( int y = yFrom; (y < dimsY) && rowsToDo; yFrom = 0, ++y, --rowsToDo )
      {
      xform->GetTransformedGridRow( dimsX, &v[0], 0, y, z );
      for ( int x = 0; x < dimsX; ++x )
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

} // namespace cmtk
