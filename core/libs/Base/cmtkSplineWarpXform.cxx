/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include "cmtkSplineWarpXform.h"

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkBitVector.h>

#include <vector>
#include <cassert>

namespace
cmtk
{

/** \addtogroup Base */
//@{

SplineWarpXform::SplineWarpXform()
{
  this->Init();
}

void SplineWarpXform::Init () 
{
  this->m_GlobalScaling = 1.0;
}

SplineWarpXform::SplineWarpXform 
( const FixedVector<3,Types::Coordinate>& domain, const Types::Coordinate delta, const AffineXform* initialXform, const bool exactDelta  )
{
  this->Init( domain, delta, initialXform, exactDelta );
}

void
SplineWarpXform
::Init
( const Self::SpaceVectorType& domain, const Types::Coordinate delta, const AffineXform* initialXform, const bool exactDelta  )
{
  this->Init();
  this->m_Domain = domain;
  if ( initialXform )
    {
    this->m_InitialAffineXform = initialXform->Clone();
    }
  else
    {
    this->m_InitialAffineXform = AffineXform::SmartPtr( NULL );
    }

  if ( exactDelta ) 
    {
    for ( int dim=0; dim<3; ++dim ) 
      {
      this->m_Spacing[dim] = delta;
      this->m_Dims[dim] = static_cast<int>( 4 + (this->m_Domain[dim] / this->m_Spacing[dim]) );
      this->m_Domain[dim] = (this->m_Dims[dim] - 3) * this->m_Spacing[dim];
      }
    } 
  else
    {
    for ( int dim=0; dim<3; ++dim )
      this->m_Dims[dim] = 2 + std::max( 2, 1+static_cast<int>( domain[dim]/delta ) );
    }
  
  this->m_NumberOfControlPoints = this->m_Dims[0] * this->m_Dims[1] * this->m_Dims[2];
  this->AllocateParameterVector( 3 * this->m_NumberOfControlPoints );
  
  this->Update( exactDelta );
  this->InitControlPoints( this->m_InitialAffineXform );
}

SplineWarpXform::SplineWarpXform
( const FixedVector<3,Types::Coordinate>& domain, const Self::ControlPointIndexType& dims, CoordinateVector::SmartPtr& parameters, const AffineXform* initialXform )
{
  this->Init();
  this->m_Domain = domain;
  this->m_Dims = dims;

  if ( initialXform )
    {
    this->m_InitialAffineXform = initialXform->Clone();
    this->m_GlobalScaling = this->m_InitialAffineXform->GetGlobalScaling();
    } 
  else
    {
    this->m_InitialAffineXform = AffineXform::SmartPtr( NULL );
    }

  this->m_NumberOfControlPoints = this->m_Dims[0] * this->m_Dims[1] * this->m_Dims[2];
  this->m_NumberOfParameters = 3 * this->m_NumberOfControlPoints;

  if ( !parameters )
    this->m_ParameterVector = CoordinateVector::SmartPtr( new CoordinateVector( this->m_NumberOfParameters ) );
  else
    this->m_ParameterVector = parameters;
  this->m_Parameters = this->m_ParameterVector->Elements;

  this->Update( false /* exactDelta */ );

  if ( !parameters )
    this->InitControlPoints( this->m_InitialAffineXform );
}

void
SplineWarpXform::InitControlPoints( const AffineXform* affineXform )
{
  Types::Coordinate *ofs = this->m_Parameters;
  Types::Coordinate pZ = -this->m_Spacing[2];
  for ( int z=0; z<this->m_Dims[2]; ++z, pZ+=this->m_Spacing[2] ) 
    {
    Types::Coordinate pY = -this->m_Spacing[1];
    for ( int y=0; y<this->m_Dims[1]; ++y, pY+=this->m_Spacing[1] ) 
      {
      Types::Coordinate pX = -this->m_Spacing[0];
      for ( int x=0; x<this->m_Dims[0]; ++x, pX+=this->m_Spacing[0], ofs+=3 ) 
	{
	ofs[0] = pX;
	ofs[1] = pY;
	ofs[2] = pZ;
	}
      }
    }
  
  if ( affineXform ) 
    {
    ofs = this->m_Parameters;
    for ( unsigned int idx = 0; idx < this->m_NumberOfControlPoints; ++idx, ofs+=3 ) 
      {
      Self::SpaceVectorType p( ofs );
      affineXform->ApplyInPlace( p );
      ofs[0] = p[0];
      ofs[1] = p[1];
      ofs[2] = p[2];
      }
    
    this->m_InverseAffineScaling = affineXform->GetScales();
    this->m_GlobalScaling = affineXform->GetGlobalScaling();
    } 
  else
    {
    this->m_InverseAffineScaling[0] = this->m_InverseAffineScaling[1] = this->m_InverseAffineScaling[2] = this->m_GlobalScaling = 1.0;
    }
}

void
SplineWarpXform::Update 
( const bool exactDelta ) 
{
  this->WarpXform::Update();

  for ( int dim=0; dim<3; ++dim ) 
    {
    assert( this->m_Dims[dim] > 3 );
    if ( exactDelta ) 
      {
      this->m_InverseSpacing[dim] = 1.0 / this->m_Spacing[dim];
      } 
    else
      {
      this->m_Spacing[dim] = this->m_Domain[dim] / (this->m_Dims[dim]-3);
      this->m_InverseSpacing[dim] = 1.0*(this->m_Dims[dim]-3) / this->m_Domain[dim];
      }
    m_Offset[dim] = -this->m_Spacing[dim];
    }
  
  int dml = 0;
  for ( int dim = 0; dim<3; ++dim )
    for ( int m = 0; m < 4; ++m )
      for ( int l = 0; l < 4; ++l, ++dml )
	GridPointOffset[dml] = dim + l * nextJ + m * nextK;
}

SplineWarpXform* 
SplineWarpXform::CloneVirtual() const
{
  SplineWarpXform *newXform = new SplineWarpXform();

  newXform->m_ParameterVector = CoordinateVector::SmartPtr( this->m_ParameterVector->Clone() );
  newXform->m_Parameters = newXform->m_ParameterVector->Elements;
  newXform->m_NumberOfParameters = this->m_NumberOfParameters;
  newXform->m_NumberOfControlPoints = this->m_NumberOfControlPoints;
  
  newXform->m_Dims = this->m_Dims;
  newXform->m_Domain = this->m_Domain;
  newXform->m_Spacing = this->m_Spacing;
  newXform->m_InverseSpacing = this->m_InverseSpacing;
  newXform->m_Offset = this->m_Offset;

  if ( this->m_ActiveFlags ) 
    {
    BitVector::SmartPtr activeFlags( this->m_ActiveFlags->Clone() );
    newXform->SetActiveFlags( activeFlags );
    }
  newXform->m_IgnoreEdge = this->m_IgnoreEdge;
  newXform->m_FastMode = this->m_FastMode;

  if ( this->m_InitialAffineXform ) 
    {
    newXform->m_InitialAffineXform = AffineXform::SmartPtr( this->m_InitialAffineXform->Clone() );
    }
  
  newXform->m_GlobalScaling = this->m_GlobalScaling;

  newXform->nextI = this->nextI;
  newXform->nextJ = this->nextJ;
  newXform->nextK = this->nextK;
  newXform->nextIJ = this->nextIJ;
  newXform->nextIK = this->nextIK;
  newXform->nextJK = this->nextJK;
  newXform->nextIJK = this->nextIJK;
  memcpy( newXform->GridPointOffset, this->GridPointOffset, sizeof( this->GridPointOffset ) );

  newXform->VolumeDims = this->VolumeDims;

  newXform->m_GridOffsets = this->m_GridOffsets;
  newXform->m_GridSpline = this->m_GridSpline;
  newXform->m_GridDerivSpline = this->m_GridDerivSpline;

  return newXform;
}

void
SplineWarpXform::Refine()
{
  if ( !this->m_ParameterVector ) return;

  Self::ControlPointIndexType newDims;
  for ( int dim=0; dim<3; ++dim ) 
    newDims[dim] = 2 * this->m_Dims[dim] - 3;

  const int newNumSamples = newDims[0] * newDims[1] * newDims[2];
  const int newNumCoefficients = 3 * newNumSamples;

  CoordinateVector::SmartPtr newParameters( new CoordinateVector( newNumCoefficients ) );
  Types::Coordinate* newCoefficients = newParameters->Elements;

  Types::Coordinate newSpacing[3];
  for ( int dim=0; dim<3; ++dim ) 
    {
    newSpacing[dim] = this->m_Domain[dim] / (newDims[dim]-3);
    }

  // no linear refinement here
  const int newNextI = 3;
  const int newNextJ = newNextI * newDims[0];
  const int newNextK = newNextJ * newDims[1];
  const int newNextIJ = newNextJ + newNextI;
  const int newNextIK = newNextK + newNextI;
  const int newNextJK = newNextK + newNextJ;
  const int newNextIJK = newNextJK + newNextI;

  Types::Coordinate level0[3][3];
  memset( level0, 0, sizeof( level0 ) );
  Types::Coordinate level1[3];
  memset( level1, 0, sizeof( level1 ) );

  Types::Coordinate *ncoeff = newCoefficients;
  for ( int z = 0; z<newDims[2]; ++z ) 
    {
    for ( int y = 0; y<newDims[1]; ++y ) 
      {
      for ( int x = 0; x<newDims[0]; ++x ) 
	{
	const int oldX = ((x+1)/2), oldY = ((y+1)/2), oldZ = ((z+1)/2);
	const int oddX = x%2, oddY = y%2, oddZ = z%2;
	
	const Types::Coordinate *coeff = m_Parameters + oldX*nextI + oldY*nextJ + oldZ*nextK;
	
	for ( int dim=0; dim<3; ++dim, ++coeff, ++ncoeff ) 
	  {	  
	  for ( int k=0; k<3; ++k ) 
	    {
	    int ofsJK = (k-1) * nextK - nextJ;
	    for ( int j=0; j<3; ++j, ofsJK += nextJ ) 
	      {
	      if ( (oddY || j) && (oddZ || k) ) 
		{
		if ( oddX ) 
		  {
		  level0[k][j] = (coeff[ofsJK-nextI] + 6 * coeff[ofsJK] + coeff[ofsJK+nextI]) / 8;
		  } 
		else
		  {
		  level0[k][j] = (coeff[ofsJK] + coeff[ofsJK+nextI]) / 2;
		  }
		}
	      }
	    }
	  
	  for ( int k=0; k<3; ++k )
	    {
	    if ( oddZ || k )
	      {
	      if ( oddY ) 
		{
		level1[k] = (level0[k][0] + 6 * level0[k][1] + level0[k][2]) / 8;
		} 
	      else
		{
		level1[k] = (level0[k][1] + level0[k][2]) / 2;
		}
	      }
	    }
	  
	  if ( oddZ ) 
	    {
	    *ncoeff = (level1[0] + 6 * level1[1] + level1[2]) / 8;
	    } 
	  else
	    {
	    *ncoeff = (level1[1] + level1[2]) / 2;
	    } 
	  }
	}
      }
    }

  this->m_NumberOfControlPoints = newNumSamples;
  this->m_NumberOfParameters = newNumCoefficients;

  this->m_ParameterVector = newParameters;
  this->m_Parameters = newCoefficients;

  this->DeleteParameterActiveFlags();
  this->m_Dims = newDims;

  for ( int dim=0; dim<3; ++dim ) 
    {
    assert( this->m_Dims[dim] > 1 );
    this->m_Spacing[dim] = newSpacing[dim];
    this->m_InverseSpacing[dim] = 1.0 / this->m_Spacing[dim];
    m_Offset[dim] = -this->m_Spacing[dim];
    }
  
  // MUST do this AFTER acutal refinement, as precomputed increments are used
  // for old grid.
  nextI = newNextI;
  nextJ = newNextJ;
  nextK = newNextK;
  nextIJ = newNextIJ;
  nextIK = newNextIK;
  nextJK = newNextJK;
  nextIJK = newNextIJK;

  int dml = 0;
  for ( int dim = 0; dim<3; ++dim )
    for ( int m = 0; m < 4; ++m )
      for ( int l = 0; l < 4; ++l, ++dml )
	GridPointOffset[dml] = dim + l * nextJ + m * nextK;

  if ( this->m_IgnoreEdge )
    this->m_IgnoreEdge = 2 * this->m_IgnoreEdge - 1;

  this->UnRegisterVolume();
}

UniformVolume::CoordinateRegionType
SplineWarpXform::GetVolumeOfInfluence( const size_t idx, const UniformVolume::CoordinateRegionType& domain, const int fastMode ) const
{
  Self::SpaceVectorType regionFrom, regionTo;
  
  int relIdx = idx / 3;

  const int xyz[3] = { ( relIdx %  this->m_Dims[0] ), 
		       ( (relIdx / this->m_Dims[0]) % this->m_Dims[1] ), 
		       ( (relIdx / this->m_Dims[0]) / this->m_Dims[1] ) };
  
  FixedVector<3,Types::Coordinate> xyzLow, xyzUp;

  if ( (fastMode==1) || (this->m_FastMode && (fastMode<0)) ) 
    {
    for ( int dim = 0; dim < 3; ++dim )
      {
      xyzLow[dim] = this->m_Spacing[dim] * std::max( 0, xyz[dim]-2 );
      xyzUp[dim] = this->m_Spacing[dim] * std::min( this->m_Dims[dim]-3, xyz[dim] );
      }
    } 
  else
    {
    for ( int dim = 0; dim < 3; ++dim )
      {
      xyzLow[dim] = this->m_Spacing[dim] * std::max( 0, xyz[dim]-3 );
      xyzUp[dim] = this->m_Spacing[dim] * std::min( this->m_Dims[dim]-2, xyz[dim]+1 );
      }
    }
  
  for ( int dim = 0; dim < 3; ++dim )
    {
    regionFrom[dim] = std::min( domain.To()[dim], std::max( xyzLow[dim], domain.From()[dim]) );
    regionTo[dim] = std::max( domain.From()[dim], std::min( xyzUp[dim], domain.To()[dim]) ); 
    }

  return UniformVolume::CoordinateRegionType( regionFrom, regionTo );
}

void
SplineWarpXform::RegisterVolumeAxis 
( const DataGrid::IndexType::ValueType dim, const Types::Coordinate delta, const Types::Coordinate origin, const int cpgDim, const Types::Coordinate invCpgSpacing,
  std::vector<int>& gIdx, std::vector<Types::Coordinate>& spline, std::vector<Types::Coordinate>& dspline )
{
  gIdx.resize( dim+1 );
  spline.resize( 4*dim );
  dspline.resize( 4*dim );

  for ( int idx=0; idx<dim; ++idx ) 
    {
    const Types::Coordinate r = invCpgSpacing * (origin + delta * idx);
    gIdx[idx] = std::min( static_cast<int>( r ), cpgDim-4 );
    const Types::Coordinate f = r - gIdx[idx];
    for ( int k = 0; k < 4; ++k ) 
      {
      spline[4*idx+k] = CubicSpline::ApproxSpline( k, f );
      dspline[4*idx+k] = CubicSpline::DerivApproxSpline( k, f );
      }
    }
  // guard element
  gIdx[dim] = -1;
}

void
SplineWarpXform::RegisterVolumePoints
( const DataGrid::IndexType& volDims, const Self::SpaceVectorType& delta, const Self::SpaceVectorType& origin )
{
  for ( int axis = 0; axis < 3; ++axis )
    this->RegisterVolumeAxis( volDims[axis], delta[axis], origin[axis], this->m_Dims[axis], this->m_InverseSpacing[axis], this->m_GridIndexes[axis], this->m_GridSpline[axis], this->m_GridDerivSpline[axis] );
  
  for ( int idx = 0; idx<volDims[0]; ++idx ) 
    this->m_GridOffsets[0][idx] = this->m_GridIndexes[0][idx] * nextI;
  for ( int idx = 0; idx<volDims[1]; ++idx ) 
    this->m_GridOffsets[1][idx] = this->m_GridIndexes[1][idx] * nextJ;
  for ( int idx = 0; idx<volDims[2]; ++idx ) 
    this->m_GridOffsets[2][idx] = this->m_GridIndexes[2][idx] * nextK;
  
  this->VolumeDims = volDims;
}

void SplineWarpXform::UnRegisterVolume()
{
  for ( int axis = 0; axis < 3; ++axis )
    {
    this->m_GridOffsets[axis].resize( 0 );
    this->m_GridSpline[axis].resize( 0 );
    this->m_GridDerivSpline[axis].resize( 0 );
    }
}

SplineWarpXform::SpaceVectorType& 
SplineWarpXform::GetDeformedControlPointPosition
( Self::SpaceVectorType& v, const int x, const int y, const int z) 
  const 
{
  // Create a pointer to the front-lower-left corner of the c.p.g. cell.
  const Types::Coordinate* coeff = m_Parameters + 3 * ( (x-1) + this->m_Dims[0] * ((y-1) + this->m_Dims[1] * (z-1)) );  
  static const Types::Coordinate spline[3] = { 1.0/6, 4.0/6, 1.0/6 };

  for ( int dim = 0; dim<3; ++dim ) 
    {
    Types::Coordinate mm = 0;
    const Types::Coordinate *coeff_mm = coeff;
    
    for ( int m = 0; m < 3; ++m ) 
      {
      Types::Coordinate ll = 0;
      const Types::Coordinate *coeff_ll = coeff_mm;
      
      // Loop over 4 c.p.g. planes in y-direction.
      for ( int l = 0; l < 3; ++l ) 
	{
	Types::Coordinate kk = 0;
	const Types::Coordinate *coeff_kk = coeff_ll;
	
	// Loop over 4 c.p.g. planes in x-direction.
	for ( int k = 0; k < 3; ++k, coeff_kk+=3 ) 
	  {
	  kk += spline[k] * (*coeff_kk);
	  }
	ll += spline[l] * kk;
	coeff_ll += nextJ;
	}	
      mm += spline[m] * ll;
      coeff_mm += nextK;
      }
    v[dim] = mm;
    ++coeff;
    }
  
  return v;
}

SplineWarpXform::SpaceVectorType
SplineWarpXform
::GetTransformedGrid 
( const int idxX, const int idxY, const int idxZ ) const
{
  Self::SpaceVectorType v;

  const Types::Coordinate* coeff = this->m_Parameters + this->m_GridOffsets[0][idxX] + this->m_GridOffsets[1][idxY] + this->m_GridOffsets[2][idxZ];
  const Types::Coordinate *spX = &this->m_GridSpline[0][idxX<<2], *spY = &this->m_GridSpline[1][idxY<<2], *spZ = &this->m_GridSpline[2][idxZ<<2];
  
  for ( int dim = 0; dim<3; ++dim ) 
    {
    Types::Coordinate mm = 0;
    const Types::Coordinate *coeff_mm = coeff;
      for ( int m = 0; m < 4; ++m )
	{
	Types::Coordinate ll = 0;
	const Types::Coordinate *coeff_ll = coeff_mm;
	for ( int l = 0; l < 4; ++l ) 
	  {
	  Types::Coordinate kk = 0;
	  const Types::Coordinate *coeff_kk = coeff_ll;
	  for ( int k = 0; k < 4; ++k, coeff_kk+=3 ) 
	    {
	    kk += spX[k] * (*coeff_kk);
	    }
	  ll += spY[l] * kk;
	  coeff_ll += nextJ;
	  }	
	mm += spZ[m] * ll;
	coeff_mm += nextK;
	}
      v[ dim ] = mm;
      ++coeff;
    }

  return v;
}

void 
SplineWarpXform::GetTransformedGridRow
( const int numPoints, Self::SpaceVectorType *const vIn, const int idxX, const int idxY, const int idxZ ) 
  const
{
  Self::SpaceVectorType *v = vIn;
  const Types::Coordinate* coeff = m_Parameters + this->m_GridOffsets[0][idxX] + this->m_GridOffsets[1][idxY] + this->m_GridOffsets[2][idxZ];
  const Types::Coordinate *spX = &this->m_GridSpline[0][idxX<<2], *spY = &this->m_GridSpline[1][idxY<<2], *spZ = &this->m_GridSpline[2][idxZ<<2];
  
  // precompute the products of B_j(v) and B_k(w) for the 4 x 4 neighborhood
  // in y- and z-direction.
  Types::Coordinate sml[16], *psml = sml;
  for ( int m = 0; m < 4; ++m )
    {
    for ( int l = 0; l < 4; ++l, ++psml )
      {
      *psml = spZ[m] * spY[l];
      }
    }
  
  // determine the number of CPG cells on our way along the row
  const int numberOfCells = (this->m_GridOffsets[0][idxX + numPoints - 1] - this->m_GridOffsets[0][idxX]) / nextI + 4;
  
  // pre-compute the contributions of all control points in y- and z-direction
  // along the way
  Types::Coordinate phiComp;
#ifdef CMTK_COMPILER_VAR_AUTO_ARRAYSIZE
  Types::Coordinate phiHat[3*numberOfCells]; // GNU compiler can have variable-sized automatic arrays
#else
  std::vector<Types::Coordinate> phiHat( 3*numberOfCells );
#endif

  const int *gpo;
  int phiIdx = 0;
  for ( int cell = 0; cell < numberOfCells; ++cell, coeff += nextI ) 
    {
    gpo = &this->GridPointOffset[0];
    for ( int dim = 0; dim < 3; ++dim, ++phiIdx ) 
      {
      phiComp = coeff[ *gpo ] * sml[0];
      ++gpo;
      for ( int ml = 1; ml < 16; ++ml, ++gpo ) 
	{
	phiComp += coeff[ *gpo ] * sml[ml];
	}
      phiHat[phiIdx] = phiComp;
      }
    }
  
  // start at the leftmost precomputed CPG cell
  int cellIdx = 0;

  // run over all points we're supposed to transform
  int i = idxX;
  for ( const int lastPoint = idxX + numPoints; i < lastPoint; ) 
    {
    // these change only when we enter a new cell
    const Types::Coordinate* phiPtr = &phiHat[3*cellIdx];
    
    // do everything inside one cell
    do 
      {
      Self::SpaceVectorType& vRef = *v;
      // compute transformed voxel by taking precomputed y- and z-contributions
      // and adding x. The loops to do this have been unrolled for increased
      // performance.
      vRef[0] = spX[0] * phiPtr[0] + spX[1] * phiPtr[3] + spX[2] * phiPtr[6] + spX[3] * phiPtr[9];
      vRef[1] = spX[0] * phiPtr[1] + spX[1] * phiPtr[4] + spX[2] * phiPtr[7] + spX[3] * phiPtr[10];
      vRef[2] = spX[0] * phiPtr[2] + spX[1] * phiPtr[5] + spX[2] * phiPtr[8] + spX[3] * phiPtr[11];
      
      // go to next voxel
      ++i;
      spX += 4;
      ++v;
      // repeat this until we leave current CPG cell.
      } 
    while ( ( this->m_GridOffsets[0][i-1] == this->m_GridOffsets[0][i] ) && ( i < lastPoint ) );
    
    // we've just left a cell -- shift index of precomputed control points
    // to the next one.
    ++cellIdx;
    }
}

Types::Coordinate
SplineWarpXform
::GetGridEnergy() const
{
  double energy = 0;

#pragma omp parallel for reduction(+:energy)
  for ( int z = 1; z<this->m_Dims[2]-1; ++z )
    {
    for ( int y = 1; y<this->m_Dims[1]-1; ++y )
      {
      for ( int x = 1; x<this->m_Dims[0]-1; ++x ) 
	{
	const Types::Coordinate* coeff = this->m_Parameters + x * nextI + y * nextJ + z * nextK;
	energy += this->GetGridEnergy( coeff );
	}
      }
    }
  
  return energy / (( this->m_Dims[0] - 2 ) * ( this->m_Dims[1] - 2 ) * ( this->m_Dims[2] - 2 ));
}

Types::Coordinate
SplineWarpXform
::GetGridEnergy( const Self::SpaceVectorType& v ) const
{
  Types::Coordinate r[3], f[3];
  int grid[3];
  
  for ( int dim = 0; dim<3; ++dim ) 
    {
    r[dim] = this->m_InverseSpacing[dim] * v[dim];
    grid[dim] = std::min( static_cast<int>( r[dim] ), this->m_Dims[dim]-4 );
    f[dim] = std::max<Types::Coordinate>( 0, std::min<Types::Coordinate>( 1.0, r[dim] - grid[dim] ) );
    }
  
  const Types::Coordinate* coeff = this->m_Parameters + 3 * ( grid[0] + this->m_Dims[0] * (grid[1] + this->m_Dims[1] * grid[2]) );

  // Matrix of one-variable second-order derivatives.
  double J[3][3];
  memset( J, 0, sizeof( J ) );

  // Matrix of mixed second-order derivatives.
  double K[3][3];
  memset( K, 0, sizeof( K ) );

  for ( int dim = 0; dim<3; ++dim ) 
    {
    const Types::Coordinate *coeff_mm = coeff;
    for ( int m = 0; m < 3; ++m ) 
      {
      Types::Coordinate llJ[3] = { 0, 0, 0 };
      Types::Coordinate llK[3] = { 0, 0, 0 };
      const Types::Coordinate *coeff_ll = coeff_mm;
      for ( int l = 0; l < 3; ++l ) 
	{
	Types::Coordinate kkJ[3] = { 0, 0, 0 };
	Types::Coordinate kkK[3] = { 0, 0, 0 };
	const Types::Coordinate *coeff_kk = coeff_ll;
	for ( int k = 0; k < 3; ++k, coeff_kk += nextI ) 
	  {
	  const Types::Coordinate tmp = CubicSpline::ApproxSpline( k, f[0] ) * (*coeff_kk);
	  kkJ[0] += CubicSpline::SecondDerivApproxSpline( k, f[0] ) * (*coeff_kk);
	  kkJ[1] += tmp;
	  kkJ[2] += tmp;
	  
	  const Types::Coordinate tmp2 = CubicSpline::DerivApproxSpline( k, f[0] ) * (*coeff_kk);
	  kkK[0] += tmp2;
	  kkK[1] += CubicSpline::ApproxSpline( k, f[0] ) * (*coeff_kk);
	  kkK[2] += tmp2;
	  }

	const Types::Coordinate tmp = CubicSpline::ApproxSpline( l, f[1] );
	llJ[0] += tmp * kkJ[0];
	llJ[1] += CubicSpline::SecondDerivApproxSpline( l, f[1] ) * kkJ[1];
	llJ[2] += tmp * kkJ[2];
	
	const Types::Coordinate tmp2 = CubicSpline::DerivApproxSpline( l, f[1] );
	llK[0] += tmp2 * kkK[0];
	llK[1] += CubicSpline::DerivApproxSpline( l, f[1] ) * kkK[1];
	llK[2] += tmp2 * kkK[2];
	coeff_ll += nextJ;
	}

      const Types::Coordinate tmp = CubicSpline::ApproxSpline( m, f[2] );
      J[0][dim] += tmp * llJ[0];
      J[1][dim] += CubicSpline::ApproxSpline( m, f[2] ) * llJ[1];
      J[2][dim] += tmp * llJ[2];
      
      const Types::Coordinate tmp2 = CubicSpline::DerivApproxSpline( m, f[2] );
      K[0][dim] += CubicSpline::ApproxSpline( m, f[2] ) * llK[0];
      K[1][dim] += tmp2 * llK[1];
      K[2][dim] += tmp2 * llK[2];
      coeff_mm += nextK;
      }
    ++coeff;
    }
  
  const double energy = 
    // single variable second-order derivatives
    MathUtil::Square( this->m_InverseSpacing[0] ) *
    ( J[0][0] * J[0][0] + J[0][1] * J[0][1] + J[0][2] * J[0][2] ) +
    MathUtil::Square( this->m_InverseSpacing[1] ) *
    ( J[1][0] * J[1][0] + J[1][1] * J[1][1] + J[1][2] * J[1][2] ) +
    MathUtil::Square( this->m_InverseSpacing[2] ) *
    ( J[2][0] * J[2][0] + J[2][1] * J[2][1] + J[2][2] * J[2][2] ) +
    // two-variable mixed derivatives
    2 * ( this->m_InverseSpacing[0] * this->m_InverseSpacing[1] *
	  ( K[0][0] * K[0][0] + K[0][1] * K[0][1] + K[0][2] * K[0][2] ) +
	  this->m_InverseSpacing[1] * this->m_InverseSpacing[2] *
	  ( K[1][0] * K[1][0] + K[1][1] * K[1][1] + K[1][2] * K[1][2] ) +
	  this->m_InverseSpacing[2] * this->m_InverseSpacing[0] *
	  ( K[2][0] * K[2][0] + K[2][1] * K[2][1] + K[2][2] * K[2][2] )
	  );
  
  return energy;
}

Types::Coordinate
SplineWarpXform
::GetGridEnergy( const Types::Coordinate *cp ) const
{
  const double   sp[3] = {  1.0/6, 2.0/3, 1.0/6 };
  const double  dsp[3] = { -1.0/2,     0, 1.0/2 };
  const double ddsp[3] = {      1,    -2,     1 };

  // Matrix of one-variable second-order derivatives.
  double J[3][3];
  memset( J, 0, sizeof( J ) );

  // Matrix of mixed second-order derivatives.
  double K[3][3];
  memset( K, 0, sizeof( K ) );

  const Types::Coordinate* coeff = cp - nextI - nextJ - nextK;
  for ( int dim = 0; dim<3; ++dim ) 
    {
    const Types::Coordinate *coeff_mm = coeff;
    for ( int m = 0; m < 3; ++m ) 
      {
      Types::Coordinate llJ[3] = { 0, 0, 0 };
      Types::Coordinate llK[3] = { 0, 0, 0 };
      const Types::Coordinate *coeff_ll = coeff_mm;
      for ( int l = 0; l < 3; ++l ) 
	{
	Types::Coordinate kkJ[3] = { 0, 0, 0 };
	Types::Coordinate kkK[3] = { 0, 0, 0 };
	const Types::Coordinate *coeff_kk = coeff_ll;
	for ( int k = 0; k < 3; ++k, coeff_kk += nextI ) 
	  {
	  const Types::Coordinate tmp = sp[k] * (*coeff_kk);
	  kkJ[0] += ddsp[k] * (*coeff_kk);
	  kkJ[1] += tmp;
	  kkJ[2] += tmp;
	  
	  const Types::Coordinate tmp2 = dsp[k] * (*coeff_kk);
	  kkK[0] += tmp2;
	  kkK[1] +=  sp[k] * (*coeff_kk);
	  kkK[2] += tmp2;
	  }
	llJ[0] +=   sp[l] * kkJ[0];
	llJ[1] += ddsp[l] * kkJ[1];
	llJ[2] +=   sp[l] * kkJ[2];
	
	llK[0] +=  dsp[l] * kkK[0];
	llK[1] +=  dsp[l] * kkK[1];
	llK[2] +=   sp[l] * kkK[2];
	coeff_ll += nextJ;
	}	
      J[0][dim] +=   sp[m] * llJ[0];
      J[1][dim] +=   sp[m] * llJ[1];
      J[2][dim] += ddsp[m] * llJ[2];
      
      K[0][dim] +=   sp[m] * llK[0];
      K[1][dim] +=  dsp[m] * llK[1];
      K[2][dim] +=  dsp[m] * llK[2];
      coeff_mm += nextK;
      }
    ++coeff;
    }
  
  const double energy = 
    // single variable second-order derivatives
    MathUtil::Square( this->m_InverseSpacing[0] ) *
    ( J[0][0] * J[0][0] + J[0][1] * J[0][1] + J[0][2] * J[0][2] ) +
    MathUtil::Square( this->m_InverseSpacing[1] ) *
    ( J[1][0] * J[1][0] + J[1][1] * J[1][1] + J[1][2] * J[1][2] ) +
    MathUtil::Square( this->m_InverseSpacing[2] ) *
    ( J[2][0] * J[2][0] + J[2][1] * J[2][1] + J[2][2] * J[2][2] ) +
    // two-variable mixed derivatives
    2 * ( this->m_InverseSpacing[0] * this->m_InverseSpacing[1] *
	  ( K[0][0] * K[0][0] + K[0][1] * K[0][1] + K[0][2] * K[0][2] ) +
	  this->m_InverseSpacing[1] * this->m_InverseSpacing[2] *
	  ( K[1][0] * K[1][0] + K[1][1] * K[1][1] + K[1][2] * K[1][2] ) +
	  this->m_InverseSpacing[2] * this->m_InverseSpacing[0] *
	  ( K[2][0] * K[2][0] + K[2][1] * K[2][1] + K[2][2] * K[2][2] )
	  );
  
  return energy;
}

void 
SplineWarpXform::GetGridEnergyDerivative
( double& lower, double& upper, const int param, const Types::Coordinate step )
  const
{
  const int controlPointIdx = param / nextI;
  const unsigned short x =  ( controlPointIdx %  this->m_Dims[0] );
  const unsigned short y = ( (controlPointIdx /  this->m_Dims[0]) % this->m_Dims[1] );
  const unsigned short z = ( (controlPointIdx /  this->m_Dims[0]) / this->m_Dims[1] );
  
  const int thisDim = param % nextI;
  const Types::Coordinate* coeff = m_Parameters + param - thisDim;
  
  double ground = 0;

  const int iFrom = std::max<int>( -1, 1-x );
  const int jFrom = std::max<int>( -1, 1-y );
  const int kFrom = std::max<int>( -1, 1-z );

  const int iTo = std::min<int>( 1, this->m_Dims[0]-2-x );
  const int jTo = std::min<int>( 1, this->m_Dims[1]-2-y );
  const int kTo = std::min<int>( 1, this->m_Dims[2]-2-z );

  for ( int k = kFrom; k < kTo; ++k )
    for ( int j = jFrom; j < jTo; ++j )
      for ( int i = iFrom; i < iTo; ++i )
	{
	ground += this->GetGridEnergy( coeff + i*nextI + j*nextJ + k*nextK );
	}

  upper = -ground;
  lower = -ground;
  
  const Types::Coordinate oldCoeff = m_Parameters[param];
  m_Parameters[param] += step;
  for ( int k = kFrom; k < kTo; ++k )
    for ( int j = jFrom; j < jTo; ++j )
      for ( int i = iFrom; i < iTo; ++i )
	upper += this->GetGridEnergy( coeff + i*nextI + j*nextJ + k*nextK );

  m_Parameters[param] = oldCoeff - step;
  for ( int k = kFrom; k < kTo; ++k )
    for ( int j = jFrom; j < jTo; ++j )
      for ( int i = iFrom; i < iTo; ++i )
	lower += this->GetGridEnergy( coeff + i*nextI + j*nextJ + k*nextK );

  m_Parameters[param] = oldCoeff;

  upper /= this->m_NumberOfControlPoints;
  lower /= this->m_NumberOfControlPoints;
}

Types::Coordinate
SplineWarpXform::GetInverseConsistencyError
( const WarpXform* inverse, const UniformVolume* volume,
  const DataGrid::RegionType* voi )
  const 
{
  Self::SpaceVectorType v, vv;
  Types::Coordinate result = 0.0;
  int count = 0;

  DataGrid::RegionType myVoi;
  const DataGrid::RegionType *pVoi = &myVoi;
  if ( voi ) 
    {
    pVoi = voi;
    } 
  else
    {
    myVoi = volume->GetWholeImageRegion();
    }

  const int dX = 1 + static_cast<int>( this->m_Spacing[0] / 2 * volume->m_Delta[AXIS_X] );
  const int dY = 1 + static_cast<int>( this->m_Spacing[1] / 2 * volume->m_Delta[AXIS_Y] );
  const int dZ = 1 + static_cast<int>( this->m_Spacing[2] / 2 * volume->m_Delta[AXIS_Z] );

  const int startX = pVoi->From()[0] - (pVoi->From()[0] % dX);
  const int startY = pVoi->From()[1] - (pVoi->From()[1] % dY);
  const int startZ = pVoi->From()[2] - (pVoi->From()[2] % dZ);

  const size_t length = pVoi->To()[0] - startX;
#ifdef CMTK_COMPILER_VAR_AUTO_ARRAYSIZE
  Self::SpaceVectorType vecArray[length];
#else
  std::vector<Self::SpaceVectorType> vecArray( length );
#endif

  for ( int z = startZ; z < pVoi->To()[2]; z += dZ ) 
    {
    for ( int y = startY; y < pVoi->To()[1]; y += dY ) 
      {
      Self::SpaceVectorType* pVec = &vecArray[0];
      this->GetTransformedGridRow( length, pVec, startX, y, z );

      for ( int x = startX; x < pVoi->To()[0]; x += dX, pVec += dX ) 
	{
	if ( inverse->InDomain( *pVec ) ) 
	  {
	  inverse->ApplyInPlace( *pVec );
	  v = volume->GetGridLocation( x, y, z );
	  v -= *pVec;
	  result += v.RootSumOfSquares();
	  ++count;
	  }
	}
      }
    }
  
  return count ? result / count : 0.0;
}

Types::Coordinate* 
SplineWarpXform::GetPureDeformation( const bool includeScale ) const
{
  const size_t numberOfParameters = this->m_NumberOfParameters;
  Types::Coordinate* points = Memory::ArrayC::Allocate<Types::Coordinate>(  numberOfParameters  );
  memcpy( points, this->m_Parameters, sizeof( *points ) * numberOfParameters );
  
  AffineXform::SmartPtr xform( this->GetInitialAffineXform()->MakeInverse() );

  if ( includeScale ) 
    {
    xform->SetScales( 1.0, 1.0, 1.0 );
    }

  Types::Coordinate* ptr = points;
  Self::SpaceVectorType u;
  for ( size_t pointIdx = 0; pointIdx < numberOfParameters / 3; ++pointIdx, ptr += 3 ) 
    {
    Self::SpaceVectorType v( ptr );
    
    // undo affine transformation component
    xform->ApplyInPlace( v );
    
    // copy the result into ouput array
    for ( unsigned int dim = 0; dim < 3; ++dim ) 
      ptr[dim] = v[dim];
    }
  
  return points;
}

}
