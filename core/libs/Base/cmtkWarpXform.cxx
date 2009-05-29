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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
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
  memcpy( Dims, dims, sizeof(Dims) );
  m_Origin.Set( 0, 0, 0 );
  
  NumberOfControlPoints = Dims[0] * Dims[1] * Dims[2];
  this->AllocateParameterVector( 3 * NumberOfControlPoints );
  this->Update();
}

void 
WarpXform::Update( const bool )
{
  nextI = 3;
  nextJ = nextI * Dims[0];
  nextK = nextJ * Dims[1];
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
WarpXform::GetGridEnergy() const
{
  return 0;
}

void
WarpXform::GetGridEnergyDerivative
( double&, double&, const int, const Types::Coordinate ) const
{
}

void 
WarpXform::GetJacobianConstraintDerivative
( double&, double&, const int, const Rect3D&, const Types::Coordinate )
  const
{
}

Types::Coordinate 
WarpXform::GetParamStep
( const size_t idx, const Types::Coordinate*, const Types::Coordinate mmStep ) const
{
  if ( !ActiveFlags.IsNull() && ! (*ActiveFlags)[idx] ) return 0;

  int controlPointIdx = idx / 3;
  unsigned short x =  ( controlPointIdx %  Dims[0] );
  unsigned short y = ( (controlPointIdx /  Dims[0]) % Dims[1] );
  unsigned short z = ( (controlPointIdx /  Dims[0]) / Dims[1] );
  
  if ( (x>=IgnoreEdge) && (x<(Dims[0]-IgnoreEdge)) && (y>=IgnoreEdge) && (y<(Dims[1]-IgnoreEdge)) && (z>=IgnoreEdge) && (z<(Dims[2]-IgnoreEdge)) ) 
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
  if ( !ActiveFlags ) 
    {
    ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  ActiveFlags->Set();
}

void
WarpXform::SetParameterActive
( const size_t index, const bool active )
{
  if ( !ActiveFlags ) 
    {
    ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  ActiveFlags->Set( index, active );
}

void 
WarpXform::SetParametersActive( const Rect3D& )
{
  if ( !ActiveFlags ) 
    {
    ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
}

void
WarpXform::SetParametersActive
( const int axis, const bool active )
{
  if ( !ActiveFlags ) 
    {
    ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  for ( unsigned int idx = (unsigned int)axis; idx < this->m_NumberOfParameters; idx += 3 )
    ActiveFlags->Set( idx, active );
}

void
WarpXform::SetParametersActive( const char* axes )
{
  if ( !ActiveFlags ) 
    {
    ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
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
  if ( !ActiveFlags ) 
    {
    ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  ActiveFlags->Reset( index );
}

int
WarpXform::GetParameterActive( const size_t index ) const
{
  if ( ActiveFlags )
    return (*ActiveFlags)[index];
  else
    return 1;
}

void
WarpXform::DeleteParameterActiveFlags()
{
  ActiveFlags = BitVector::SmartPtr::Null;
}

void
WarpXform::Regularize( const int weight0, const int weight1 )
{
  Types::Coordinate* ncoeff = Memory::AllocateArray<Types::Coordinate>( this->m_NumberOfParameters );
  Types::Coordinate* pcoeff = ncoeff;

  Types::Coordinate* coeff = this->m_Parameters;
  for ( int z = 0; z<Dims[2]; ++z )
    for ( int y = 0; y<Dims[1]; ++y )
      for ( int x = 0; x<Dims[0]; ++x ) 
	{
	for ( int dim = 0; dim<3; ++dim, ++coeff, ++pcoeff ) 
	  {	  
	  *pcoeff = 0;
	  int count = 0;
	  for ( int k=-1; k<2; ++k )
	    for ( int j=-1; j<2; ++j )
	      for ( int i=-1; i<2; ++i )
		if ( ((x+i)>=0) && ((y+j)>=0) && ((z+k)>=0) && ((x+i)<Dims[0]) && ((y+j)<Dims[1]) && ((z+k)<Dims[2]) ) 
		  {
		  int weight = (i?weight1:weight0) * (j?weight1:weight0) * (k?weight1:weight0);
		  *pcoeff += weight * coeff[ i*nextI + j*nextJ + k*nextK ];
		  count += weight;
		  }
	  *pcoeff /= count;
	  }
	}
  
  memcpy( coeff, ncoeff, sizeof( *coeff ) * this->m_NumberOfParameters );
  delete[] ncoeff;
}

void
WarpXform::RegisterVolume( const UniformVolume *volume )
{
  this->RegisterVolumePoints( volume->Dims, volume->Delta, volume->m_Origin.XYZ );
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
  if ( InitialAffineXform ) 
    {
    change = *(InitialAffineXform->GetInverse());
    }

  // second, append current initial affine transformation
  if ( newAffineXform )
    change.Concat( *newAffineXform );

  // apply effective change to all control points.
  Types::Coordinate *coeff = this->m_Parameters;
  for ( unsigned int idx = 0; idx < NumberOfControlPoints; ++idx, coeff+=3 ) 
    {
    Vector3D p( coeff );
    change.ApplyInPlaceNonVirtual( p );
    coeff[0] = p[0];
    coeff[1] = p[1];
    coeff[2] = p[2];
    }

  // Finally, copy new transformation. We want to create a new object here
  // if the current transformation is linked somewhere else.
  if ( newAffineXform )
    {
    InitialAffineXform = AffineXform::SmartPtr( newAffineXform->Clone() );
    }
  else
    {
    InitialAffineXform = AffineXform::SmartPtr( new AffineXform );
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
    affineXform->ApplyInPlaceNonVirtual( p );
    coeff[0] = p[0];
    coeff[1] = p[1];
    coeff[2] = p[2];
    }

  // Finally, generate combined affine transformation. We want to create a new
  // object here if the current transformation is linked somewhere else.
  if ( InitialAffineXform.GetReferenceCount() != 1 )
    InitialAffineXform = AffineXform::SmartPtr( InitialAffineXform->Clone() );
  InitialAffineXform->Concat( *affineXform );
}

} // namespace cmtk
