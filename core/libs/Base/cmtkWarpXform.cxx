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

#include <cmtkWarpXform.h>

#include <cmtkMathUtil.h>
#include <cmtkUniformVolume.h>

#include <string.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
WarpXform::InitGrid
( const Types::Coordinate domain[3], const int dims[3] )
{
  memcpy( Domain, domain, sizeof(Domain) );
  memcpy( this->m_Dims, dims, sizeof(this->m_Dims) );
  m_Offset.Set( 0, 0, 0 );
  
  NumberOfControlPoints = this->m_Dims[0] * this->m_Dims[1] * this->m_Dims[2];
  this->AllocateParameterVector( 3 * NumberOfControlPoints );
  this->Update();
}

void 
WarpXform::Update( const bool )
{
  nextI = 3;
  nextJ = nextI * this->m_Dims[0];
  nextK = nextJ * this->m_Dims[1];
  nextIJ = nextJ + nextI;
  nextIK = nextK + nextI;
  nextJK = nextK + nextJ;
  nextIJK = nextJK + nextI;
}

void
WarpXform::GetDerivativeLandmarksMSD
( double& lowerMSD, double& upperMSD, const MatchedLandmarkList* ll,
  const unsigned int idx, const Types::Coordinate step )
{
  upperMSD = lowerMSD = 0;

  Types::Coordinate pOld = this->m_Parameters[idx];

  this->m_Parameters[idx] += step;
  MatchedLandmarkList::const_iterator it = ll->begin();
  while ( it != ll->end() ) 
    {
    Vector3D source( (*it)->GetLocation() );
    Vector3D target( (*it)->GetTargetLocation() );
    this->ApplyInPlace( source );
    upperMSD += (source - target).Square();
    ++it;
    }
  
  this->m_Parameters[idx] = pOld - step;
  it = ll->begin();
  while ( it != ll->end() ) 
    {
    Vector3D source( (*it)->GetLocation() );
    Vector3D target( (*it)->GetTargetLocation() );
    this->ApplyInPlace( source );
    lowerMSD += (source - target).Square();
    ++it;
    }
  this->m_Parameters[idx] = pOld;
  
  upperMSD /= ll->size();
  lowerMSD /= ll->size();
}

Types::Coordinate
WarpXform::GetInverseConsistencyError
( const WarpXform* inverse, const UniformVolume* volume, const Rect3D* voi ) const 
{
  Vector3D v, vv;
  Types::Coordinate result = 0.0;
  int count = 0;

  Rect3D myVoi;
  const Rect3D *pVoi = &myVoi;
  if ( voi ) 
    {
    pVoi = voi;
    } 
  else
    {
    myVoi.startX = myVoi.startY = myVoi.startZ = 0;
    myVoi.endX = volume->GetDims( AXIS_X );
    myVoi.endY = volume->GetDims( AXIS_Y );
    myVoi.endZ = volume->GetDims( AXIS_Z );
    }

  for ( int z = pVoi->startZ; z < pVoi->endZ; ++z )
    for ( int y = pVoi->startY; y < pVoi->endY; ++y )
      for ( int x = pVoi->startX; x < pVoi->endX; ++x ) 
	{
	volume->GetGridLocation( v, x, y, z );
	vv = v;
	this->ApplyInPlace( vv );
	if ( inverse->InDomain( vv ) ) 
	  {
	  inverse->ApplyInPlace( vv );
	  v -= vv;
	  result += v.EuclidNorm();
	  ++count;
	  }
	}
  
  return count ? result / count : 0.0;
}

void
WarpXform::GetDerivativeInverseConsistencyError
( double& lower, double& upper, const WarpXform* inverse,
  const UniformVolume* volume, const Rect3D* voi, 
  const unsigned int idx, const Types::Coordinate step )
{
  const Types::Coordinate pOld = this->m_Parameters[idx];
  upper = lower = (-this->GetInverseConsistencyError( inverse, volume, voi ));

  this->m_Parameters[idx] += step;
  upper += this->GetInverseConsistencyError( inverse, volume, voi );

  this->m_Parameters[idx] = pOld - step;
  lower+= this->GetInverseConsistencyError( inverse, volume, voi );

  this->m_Parameters[idx] = pOld;
}

Types::Coordinate 
WarpXform::GetParamStep
( const size_t idx, const Types::Coordinate*, const Types::Coordinate mmStep ) const
{
  if ( this->m_ActiveFlags && ! (*this->m_ActiveFlags)[idx] ) return 0;

  int controlPointIdx = idx / 3;
  unsigned short x =  ( controlPointIdx %  this->m_Dims[0] );
  unsigned short y = ( (controlPointIdx /  this->m_Dims[0]) % this->m_Dims[1] );
  unsigned short z = ( (controlPointIdx /  this->m_Dims[0]) / this->m_Dims[1] );
  
  if ( (x>=this->m_IgnoreEdge) && (x<(this->m_Dims[0]-this->m_IgnoreEdge)) && 
       (y>=this->m_IgnoreEdge) && (y<(this->m_Dims[1]-this->m_IgnoreEdge)) && 
       (z>=this->m_IgnoreEdge) && (z<(this->m_Dims[2]-this->m_IgnoreEdge)) ) 
    {
    return mmStep;
    } 
  else
    {
    return 0;
    }
}

void 
WarpXform::SetParameterActive()
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  this->m_ActiveFlags->Set();
}

void
WarpXform::SetParameterActive
( const size_t index, const bool active )
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  this->m_ActiveFlags->Set( index, active );
}

void 
WarpXform::SetParametersActive( const Rect3D& )
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
}

void
WarpXform::SetParametersActive
( const int axis, const bool active )
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  for ( unsigned int idx = (unsigned int)axis; idx < this->m_NumberOfParameters; idx += 3 )
    this->m_ActiveFlags->Set( idx, active );
}

void
WarpXform::SetParametersActive( const char* axes )
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  if ( axes ) 
    {
    if ( strchr( axes, 'x' ) || strchr( axes, 'X' ) )
      this->SetParametersActive( AXIS_X );
    if ( strchr( axes, 'y' ) || strchr( axes, 'Y' ) )
      this->SetParametersActive( AXIS_Y );
    if ( strchr( axes, 'z' ) || strchr( axes, 'Z' ) )
      this->SetParametersActive( AXIS_Z );
    }
}

void
WarpXform::SetParameterInactive( const size_t index )
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  this->m_ActiveFlags->Reset( index );
}

int
WarpXform::GetParameterActive( const size_t index ) const
{
  if ( this->m_ActiveFlags )
    return (*this->m_ActiveFlags)[index];
  else
    return 1;
}

void
WarpXform::DeleteParameterActiveFlags()
{
  this->m_ActiveFlags = BitVector::SmartPtr::Null;
}

void
WarpXform::RegisterVolume( const UniformVolume *volume )
{
  this->RegisterVolumePoints( volume->m_Dims, volume->m_Delta, volume->m_Offset.XYZ );
}

void 
WarpXform::SetIncompressibilityMap
( DataGrid::SmartPtr& incompressibility )
{
  IncompressibilityMap = incompressibility;
}

void
WarpXform::ReplaceInitialAffine( const AffineXform* newAffineXform )
{
  AffineXform change;

  // first, get inverse of current initial affine transformation
  if ( this->m_InitialAffineXform ) 
    {
    change = *(this->m_InitialAffineXform->GetInverse());
    }

  // second, append current initial affine transformation
  if ( newAffineXform )
    change.Concat( *newAffineXform );

  // apply effective change to all control points.
  Types::Coordinate *coeff = this->m_Parameters;
  for ( unsigned int idx = 0; idx < NumberOfControlPoints; ++idx, coeff+=3 ) 
    {
    Vector3D p( coeff );
    change.ApplyInPlace( p );
    coeff[0] = p[0];
    coeff[1] = p[1];
    coeff[2] = p[2];
    }

  // Finally, copy new transformation. We want to create a new object here
  // if the current transformation is linked somewhere else.
  if ( newAffineXform )
    {
    this->m_InitialAffineXform = AffineXform::SmartPtr( newAffineXform->Clone() );
    }
  else
    {
    this->m_InitialAffineXform = AffineXform::SmartPtr( new AffineXform );
    }
}

void
WarpXform::ConcatAffine( const AffineXform* affineXform )
{
  // apply effective change to all control points.
  Types::Coordinate *coeff = this->m_Parameters;
  for ( unsigned int idx = 0; idx < NumberOfControlPoints; ++idx, coeff+=3 ) 
    {
    Vector3D p( coeff );
    affineXform->ApplyInPlace( p );
    coeff[0] = p[0];
    coeff[1] = p[1];
    coeff[2] = p[2];
    }

  // Finally, generate combined affine transformation. We want to create a new
  // object here if the current transformation is linked somewhere else.
  if ( this->m_InitialAffineXform.GetReferenceCount() != 1 )
    this->m_InitialAffineXform = AffineXform::SmartPtr( this->m_InitialAffineXform->Clone() );
  this->m_InitialAffineXform->Concat( *affineXform );
}

} // namespace cmtk
