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
#include <System/cmtkThreadPool.h>

#include <Base/cmtkRegionIndexIterator.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
SplineWarpXform::GetJacobianRow
( CoordinateMatrix3x3 *const array, const int x, const int y, const int z, const size_t numberOfPoints ) 
  const
{
  const Types::Coordinate *spX = &this->m_GridSpline[0][x<<2], *spY = &this->m_GridSpline[1][y<<2], *spZ = &this->m_GridSpline[2][z<<2];
  const Types::Coordinate *dspX = &this->m_GridDerivSpline[0][x<<2], *dspY = &this->m_GridDerivSpline[1][y<<2], *dspZ = &this->m_GridDerivSpline[2][z<<2];
  const Types::Coordinate *coeff = this->m_Parameters + this->m_GridOffsets[0][x] + this->m_GridOffsets[1][y] + this->m_GridOffsets[2][z];

  // precompute the products of B_j(v) and B_k(w) for the 4 x 4 neighborhood
  // in y- and z-direction.
  Types::Coordinate smlX[16], smlY[16], smlZ[16];
  for ( int i = 0, m = 0; m < 4; ++m )
    {
    for ( int l = 0; l < 4; ++l, ++i )
      {
      smlX[i] =  spZ[m] *  spY[l];
      smlY[i] =  spZ[m] * dspY[l];
      smlZ[i] = dspZ[m] *  spY[l];
      }
    }
  
  // determine the number of CPG cells on our way along the row
  const int numberOfCells = (this->m_GridOffsets[0][x + numberOfPoints - 1] - this->m_GridOffsets[0][x]) / nextI + 4;
  
  // pre-compute the contributions of all control points in y- and z-direction
  // along the way
  Types::Coordinate phiCompX, phiCompY, phiCompZ;

#ifdef CMTK_COMPILER_VAR_AUTO_ARRAYSIZE
  Types::Coordinate phiHatX[3*numberOfCells];
  Types::Coordinate phiHatY[3*numberOfCells];
  Types::Coordinate phiHatZ[3*numberOfCells];
#else
  std::vector<Types::Coordinate> phiHatX(3*numberOfCells);
  std::vector<Types::Coordinate> phiHatY(3*numberOfCells);
  std::vector<Types::Coordinate> phiHatZ(3*numberOfCells);
#endif

  const int *gpo;
  int phiIdx = 0;
  for ( int cell = 0; cell < numberOfCells; ++cell, coeff += nextI ) 
    {
    gpo = &GridPointOffset[0];
    for ( int dim = 0; dim < 3; ++dim, ++phiIdx ) 
      {
      phiCompX = phiCompY = phiCompZ = 0;
      for ( int ml = 0; ml < 16; ++ml, ++gpo ) 
	{
	phiCompX += coeff[ *gpo ] * smlX[ml];
	phiCompY += coeff[ *gpo ] * smlY[ml];
	phiCompZ += coeff[ *gpo ] * smlZ[ml];
	}
      phiHatX[phiIdx] = phiCompX;
      phiHatY[phiIdx] = phiCompY;
      phiHatZ[phiIdx] = phiCompZ;
      }
    }
  
  // start at the leftmost precomputed CPG cell
  int cellIdx = 0;

  CoordinateMatrix3x3 J;
  // run over all points we're supposed to transform
  int i = x;
  for ( const int lastPoint = x + numberOfPoints; i < lastPoint; ) 
    {
    // these change only when we enter a new cell
#ifdef CMTK_COMPILER_VAR_AUTO_ARRAYSIZE
    const Types::Coordinate* phiPtrX = phiHatX + 3*cellIdx;
    const Types::Coordinate* phiPtrY = phiHatY + 3*cellIdx;
    const Types::Coordinate* phiPtrZ = phiHatZ + 3*cellIdx;
#else
    const Types::Coordinate* phiPtrX = &phiHatX[3*cellIdx];
    const Types::Coordinate* phiPtrY = &phiHatY[3*cellIdx];
    const Types::Coordinate* phiPtrZ = &phiHatZ[3*cellIdx];
#endif

    // do everything inside one cell
    do 
      {
      // compute transformed voxel by taking precomputed y- and z-contributions
      // and adding x. The loops to do this have been unrolled for increased
      // performance.
      J[0][0] = this->m_InverseSpacing[0] * ( dspX[0] * phiPtrX[0] + dspX[1] * phiPtrX[3] + dspX[2] * phiPtrX[6] + dspX[3] * phiPtrX[9] );
      J[0][1] = this->m_InverseSpacing[0] * ( dspX[0] * phiPtrX[1] + dspX[1] * phiPtrX[4] + dspX[2] * phiPtrX[7] + dspX[3] * phiPtrX[10] );
      J[0][2] = this->m_InverseSpacing[0] * ( dspX[0] * phiPtrX[2] + dspX[1] * phiPtrX[5] + dspX[2] * phiPtrX[8] + dspX[3] * phiPtrX[11] );

      J[1][0] = this->m_InverseSpacing[1] * ( spX[0] * phiPtrY[0] + spX[1] * phiPtrY[3] + spX[2] * phiPtrY[6] + spX[3] * phiPtrY[9] );
      J[1][1] = this->m_InverseSpacing[1] * ( spX[0] * phiPtrY[1] + spX[1] * phiPtrY[4] + spX[2] * phiPtrY[7] + spX[3] * phiPtrY[10] );
      J[1][2] = this->m_InverseSpacing[1] * ( spX[0] * phiPtrY[2] + spX[1] * phiPtrY[5] + spX[2] * phiPtrY[8] + spX[3] * phiPtrY[11] );

      J[2][0] = this->m_InverseSpacing[2] * ( spX[0] * phiPtrZ[0] + spX[1] * phiPtrZ[3] + spX[2] * phiPtrZ[6] + spX[3] * phiPtrZ[9] );
      J[2][1] = this->m_InverseSpacing[2] * ( spX[0] * phiPtrZ[1] + spX[1] * phiPtrZ[4] + spX[2] * phiPtrZ[7] + spX[3] * phiPtrZ[10] );
      J[2][2] = this->m_InverseSpacing[2] * ( spX[0] * phiPtrZ[2] + spX[1] * phiPtrZ[5] + spX[2] * phiPtrZ[8] + spX[3] * phiPtrZ[11] );

      array[i-x].Set( &J[0][0] );

      // go to next voxel
      ++i;
      spX += 4;
      dspX += 4;
      // repeat this until we leave current CPG cell.
    } while ( ( this->m_GridOffsets[0][i-1] == this->m_GridOffsets[0][i] ) && ( i < lastPoint ) );

    // we've just left a cell -- shift index of precomputed control points
    // to the next one.
    ++cellIdx;
  }
}

CoordinateMatrix3x3
SplineWarpXform::GetJacobianAtControlPoint
( const Types::Coordinate* cp ) const
{
  CoordinateMatrix3x3 J = CoordinateMatrix3x3::Zero();
  
  const double  sp[3] = {  1.0/6, 2.0/3, 1.0/6 };
  const double dsp[3] = { -1.0/2,     0, 1.0/2 };

  const Types::Coordinate* coeff = cp - nextI - nextJ - nextK;
  for ( int dim = 0; dim<3; ++dim ) {
    const Types::Coordinate *coeff_mm = coeff;
    for ( int m = 0; m < 3; ++m ) {
      Types::Coordinate ll[3] = { 0, 0, 0 };
      const Types::Coordinate *coeff_ll = coeff_mm;
      for ( int l = 0; l < 3; ++l ) {
	Types::Coordinate kk[3] = { 0, 0, 0 };
	const Types::Coordinate *coeff_kk = coeff_ll;
	for ( int k = 0; k < 3; ++k, coeff_kk += nextI ) {
	  kk[0] += dsp[k] * (*coeff_kk);
	  kk[1] +=  sp[k] * (*coeff_kk);
	  kk[2] +=  sp[k] * (*coeff_kk);
	}
	ll[0] +=  sp[l] * kk[0];
	ll[1] += dsp[l] * kk[1];
	ll[2] +=  sp[l] * kk[2];
	coeff_ll += nextJ;
      }	
      J[0][dim] +=  sp[m] * ll[0];
      J[1][dim] +=  sp[m] * ll[1];
      J[2][dim] += dsp[m] * ll[2];
      coeff_mm += nextK;
    }
    ++coeff;
  }
  for ( int i = 0; i<3; ++i ) 
    {
    for ( int j = 0; j<3; ++j )
      J[i][j] *= this->m_InverseSpacing[i];
    }

  return J;
}

CoordinateMatrix3x3
SplineWarpXform::GetJacobian
( const Self::SpaceVectorType& v ) const
{
  Types::Coordinate r[3], f[3];
  int grid[3];
  
  for ( int dim = 0; dim<3; ++dim ) 
    {
    r[dim] = this->m_InverseSpacing[dim] * v[dim];
    grid[dim] = std::min( static_cast<int>( r[dim] ), this->m_Dims[dim]-4 );
    f[dim] = std::max<Types::Coordinate>( 0, std::min<Types::Coordinate>( 1.0, r[dim] - grid[dim] ) );
    }

  CoordinateMatrix3x3 J = CoordinateMatrix3x3::Zero();
  const Types::Coordinate* coeff = this->m_Parameters + 3 * ( grid[0] + this->m_Dims[0] * (grid[1] + this->m_Dims[1] * grid[2]) );
  
  // loop over the three components of the coordinate transformation function,
  // x, y, z.
  for ( int dim = 0; dim<3; ++dim, ++coeff ) 
    {
    const Types::Coordinate *coeff_mm = coeff;
    for ( int m = 0; m < 4; ++m, coeff_mm += nextK ) 
      {
      Types::Coordinate ll[3] = { 0, 0, 0 };
      const Types::Coordinate *coeff_ll = coeff_mm;
      for ( int l = 0; l < 4; ++l, coeff_ll += nextJ ) 
	{
	Types::Coordinate kk[3] = { 0, 0, 0 };
	const Types::Coordinate *coeff_kk = coeff_ll;
	
	for ( int k = 0; k < 4; ++k, coeff_kk+=3 ) 
	  {
	  // dT / dx
	  kk[0] += CubicSpline::DerivApproxSpline( k, f[0] ) * (*coeff_kk);
	  // dT / dy
	  const Types::Coordinate tmp = CubicSpline::ApproxSpline( k, f[0] ) * (*coeff_kk);
	  kk[1] += tmp;
	  // dT / dz
	  kk[2] += tmp;
	  }
	
	// dT / dx
	const Types::Coordinate tmp = CubicSpline::ApproxSpline( l, f[1] );
	ll[0] += tmp * kk[0];
	// dT / dy
	ll[1] += CubicSpline::DerivApproxSpline( l, f[1] ) * kk[1];
	// dT / dz
	ll[2] += tmp * kk[2];
	}

      // dT / dx
      const Types::Coordinate tmp = CubicSpline::ApproxSpline( m, f[2] );
      J[dim][0] += tmp * ll[0];
      // dT / dy
      J[dim][1] += tmp * ll[1];
      // dT / dz
      J[dim][2] += CubicSpline::DerivApproxSpline( m, f[2] ) * ll[2];
      }
    }

  // scale with grid spacing to normalize Jacobian (chain rule of derivation)
  for ( int i = 0; i<3; ++i ) 
    {
    for ( int j = 0; j<3; ++j )
      J[i][j] *= this->m_InverseSpacing[i];
    }

  return J;
}

Types::Coordinate
SplineWarpXform::GetJacobianDeterminant
( const int x, const int y, const int z ) const
{
  const Types::Coordinate *spX = &this->m_GridSpline[0][x<<2], *spY = &this->m_GridSpline[1][y<<2], *spZ = &this->m_GridSpline[2][z<<2];
  const Types::Coordinate *dspX = &this->m_GridDerivSpline[0][x<<2], *dspY = &this->m_GridDerivSpline[1][y<<2], *dspZ = &this->m_GridDerivSpline[2][z<<2];
  const Types::Coordinate *coeff = this->m_Parameters + this->m_GridOffsets[0][x] + this->m_GridOffsets[1][y] + this->m_GridOffsets[2][z];
  
  double J[3][3];
  memset( J, 0, sizeof( J ) );
  for ( int dim = 0; dim<3; ++dim ) 
    {
    const Types::Coordinate *coeff_mm = coeff;
    for ( int m = 0; m < 4; ++m ) 
      {
      Types::Coordinate ll[3] = { 0, 0, 0 };
      const Types::Coordinate *coeff_ll = coeff_mm;
      for ( int l = 0; l < 4; ++l ) 
	{
	Types::Coordinate kk[3] = { 0, 0, 0 };
	const Types::Coordinate *coeff_kk = coeff_ll;
	for ( int k = 0; k < 4; ++k, coeff_kk+=3 ) 
	  {
	  kk[0] += dspX[k] * (*coeff_kk);
	  kk[1] += spX[k];
	  kk[2] += spX[k];
	  }
	ll[0] += spY[l] * kk[0];
	ll[1] += dspY[l] * kk[1];
	ll[2] += spY[l] * kk[2];
	coeff_ll += nextJ;
	}	
      J[0][dim] += spZ[m] * ll[0];
      J[1][dim] += spZ[m] * ll[1];
      J[2][dim] += dspZ[m] * ll[2];
      coeff_mm += nextK;
      }
    ++coeff;
    }
  
  return this->m_InverseSpacing[0] * this->m_InverseSpacing[1] * this->m_InverseSpacing[2] * 
    ( J[0][0] * (J[1][1]*J[2][2] - J[1][2]*J[2][1]) + 
      J[0][1] * (J[1][2]*J[2][0] - J[1][0]*J[2][2]) + 
      J[0][2] * (J[1][0]*J[2][1] - J[1][1]*J[2][0]) );
}

Types::Coordinate
SplineWarpXform::GetJacobianDeterminant
( const Self::SpaceVectorType& v ) const
{
  Types::Coordinate r[3], f[3];
  int grid[3];
  
  double J[3][3];
  memset( J, 0, sizeof( J ) );

  for ( int dim = 0; dim<3; ++dim ) 
    {
    r[dim] = this->m_InverseSpacing[dim] * v[dim];
    grid[dim] = std::min( static_cast<int>( r[dim] ), this->m_Dims[dim]-4 );
    f[dim] = std::max<Types::Coordinate>( 0, std::min<Types::Coordinate>( 1.0, r[dim] - grid[dim] ) );
    }
  
  const Types::Coordinate* coeff = this->m_Parameters + 3 * ( grid[0] + this->m_Dims[0] * (grid[1] + this->m_Dims[1] * grid[2]) );
  
  for ( int dim = 0; dim<3; ++dim ) 
    {
    const Types::Coordinate *coeff_mm = coeff;
    for ( int m = 0; m < 4; ++m ) 
      {
      Types::Coordinate ll[3] = { 0, 0, 0 };
      const Types::Coordinate *coeff_ll = coeff_mm;
      for ( int l = 0; l < 4; ++l ) 
	{
	Types::Coordinate kk[3] = { 0, 0, 0 };
	const Types::Coordinate *coeff_kk = coeff_ll;
	for ( int k = 0; k < 4; ++k, coeff_kk+=3 ) 
	  {
	  kk[0] += CubicSpline::DerivApproxSpline( k, f[0] ) * (*coeff_kk);
	  const Types::Coordinate tmp = CubicSpline::ApproxSpline( k, f[0] ) * (*coeff_kk);
	  kk[1] += tmp;
	  kk[2] += tmp;
	  }
	const Types::Coordinate tmp = CubicSpline::ApproxSpline( l, f[1] );
	ll[0] += tmp * kk[0];
	ll[1] += CubicSpline::DerivApproxSpline( l, f[1] ) * kk[1];
	ll[2] += tmp * kk[2];
	coeff_ll += nextJ;
	}	
      const Types::Coordinate tmp = CubicSpline::ApproxSpline( m, f[2] );
      J[0][dim] += tmp * ll[0];
      J[1][dim] += tmp * ll[1];
      J[2][dim] += CubicSpline::DerivApproxSpline( m, f[2] ) * ll[2];
      coeff_mm += nextK;
      }
    ++coeff;
    }
  
  return this->m_InverseSpacing[0] * this->m_InverseSpacing[1] * this->m_InverseSpacing[2] * 
    ( J[0][0] * (J[1][1]*J[2][2] - J[1][2]*J[2][1]) + 
      J[0][1] * (J[1][2]*J[2][0] - J[1][0]*J[2][2]) + 
      J[0][2] * (J[1][0]*J[2][1] - J[1][1]*J[2][0]) );
}

void
SplineWarpXform::GetJacobianDeterminantRow
( double *const values, const int x, const int y, const int z, 
  const size_t numberOfPoints ) 
  const
{
  const Types::Coordinate *spX = &this->m_GridSpline[0][x<<2], *spY = &this->m_GridSpline[1][y<<2], *spZ = &this->m_GridSpline[2][z<<2];
  const Types::Coordinate *dspX = &this->m_GridDerivSpline[0][x<<2], *dspY = &this->m_GridDerivSpline[1][y<<2], *dspZ = &this->m_GridDerivSpline[2][z<<2];
  const Types::Coordinate *coeff = this->m_Parameters + this->m_GridOffsets[0][x] + this->m_GridOffsets[1][y] + this->m_GridOffsets[2][z];

  const Types::Coordinate globalInverseSpacing = this->m_InverseSpacing[0] * this->m_InverseSpacing[1] * this->m_InverseSpacing[2];

  // precompute the products of B_j(v) and B_k(w) for the 4 x 4 neighborhood
  // in y- and z-direction.
  Types::Coordinate smlX[16], smlY[16], smlZ[16];
  for ( int i = 0, m = 0; m < 4; ++m )
    {
    for ( int l = 0; l < 4; ++l, ++i )
      {
      smlX[i] =  spZ[m] *  spY[l];
      smlY[i] =  spZ[m] * dspY[l];
      smlZ[i] = dspZ[m] *  spY[l];
      }
    }
  
  // determine the number of CPG cells on our way along the row
  const int numberOfCells = (this->m_GridOffsets[0][x + numberOfPoints - 1] - this->m_GridOffsets[0][x]) / nextI + 4;
  
  // pre-compute the contributions of all control points in y- and z-direction
  // along the way
  Types::Coordinate phiCompX, phiCompY, phiCompZ;

#ifdef CMTK_COMPILER_VAR_AUTO_ARRAYSIZE
  Types::Coordinate phiHatX[3*numberOfCells];
  Types::Coordinate phiHatY[3*numberOfCells];
  Types::Coordinate phiHatZ[3*numberOfCells];
#else
  std::vector<Types::Coordinate> phiHatX(3*numberOfCells);
  std::vector<Types::Coordinate> phiHatY(3*numberOfCells);
  std::vector<Types::Coordinate> phiHatZ(3*numberOfCells);
#endif

  const int *gpo;
  int phiIdx = 0;
  for ( int cell = 0; cell < numberOfCells; ++cell, coeff += nextI ) 
    {
    gpo = &GridPointOffset[0];
    for ( int dim = 0; dim < 3; ++dim, ++phiIdx ) 
      {
      phiCompX = phiCompY = phiCompZ = 0;
      for ( int ml = 0; ml < 16; ++ml, ++gpo ) 
	{
	phiCompX += coeff[ *gpo ] * smlX[ml];
	phiCompY += coeff[ *gpo ] * smlY[ml];
	phiCompZ += coeff[ *gpo ] * smlZ[ml];
	}
      phiHatX[phiIdx] = phiCompX;
      phiHatY[phiIdx] = phiCompY;
      phiHatZ[phiIdx] = phiCompZ;
      }
    }
  
  // start at the leftmost precomputed CPG cell
  int cellIdx = 0;

  Types::Coordinate JXX, JXY, JXZ, JYX, JYY, JYZ, JZX, JZY, JZZ;
  // run over all points we're supposed to transform
  int i = x;
  for ( const int lastPoint = x + numberOfPoints; i < lastPoint; ) 
    {
    // these change only when we enter a new cell
#ifdef CMTK_COMPILER_VAR_AUTO_ARRAYSIZE
    const Types::Coordinate* phiPtrX = phiHatX + 3*cellIdx;
    const Types::Coordinate* phiPtrY = phiHatY + 3*cellIdx;
    const Types::Coordinate* phiPtrZ = phiHatZ + 3*cellIdx;
#else
    const Types::Coordinate* phiPtrX = &phiHatX[3*cellIdx];
    const Types::Coordinate* phiPtrY = &phiHatY[3*cellIdx];
    const Types::Coordinate* phiPtrZ = &phiHatZ[3*cellIdx];
#endif

    // do everything inside one cell
    do 
      {
      // compute transformed voxel by taking precomputed y- and z-contributions
      // and adding x. The loops to do this have been unrolled for increased
      // performance.
      JXX = dspX[0] * phiPtrX[0] + dspX[1] * phiPtrX[3] + dspX[2] * phiPtrX[6] + dspX[3] * phiPtrX[9];
      JXY = dspX[0] * phiPtrX[1] + dspX[1] * phiPtrX[4] + dspX[2] * phiPtrX[7] + dspX[3] * phiPtrX[10];
      JXZ = dspX[0] * phiPtrX[2] + dspX[1] * phiPtrX[5] + dspX[2] * phiPtrX[8] + dspX[3] * phiPtrX[11];

      JYX = spX[0] * phiPtrY[0] + spX[1] * phiPtrY[3] + spX[2] * phiPtrY[6] + spX[3] * phiPtrY[9];
      JYY = spX[0] * phiPtrY[1] + spX[1] * phiPtrY[4] + spX[2] * phiPtrY[7] + spX[3] * phiPtrY[10];
      JYZ = spX[0] * phiPtrY[2] + spX[1] * phiPtrY[5] + spX[2] * phiPtrY[8] + spX[3] * phiPtrY[11];

      JZX = spX[0] * phiPtrZ[0] + spX[1] * phiPtrZ[3] + spX[2] * phiPtrZ[6] + spX[3] * phiPtrZ[9];
      JZY = spX[0] * phiPtrZ[1] + spX[1] * phiPtrZ[4] + spX[2] * phiPtrZ[7] + spX[3] * phiPtrZ[10];
      JZZ = spX[0] * phiPtrZ[2] + spX[1] * phiPtrZ[5] + spX[2] * phiPtrZ[8] + spX[3] * phiPtrZ[11];

      values[i-x] = globalInverseSpacing * 
	( JXX * (JYY*JZZ - JYZ*JZY) + JXY * (JYZ*JZX - JYX*JZZ) + JXZ * (JYX*JZY - JYY*JZX) );

      // go to next voxel
      ++i;
      spX += 4;
      dspX += 4;
      // repeat this until we leave current CPG cell.
    } while ( ( this->m_GridOffsets[0][i-1] == this->m_GridOffsets[0][i] ) && ( i < lastPoint ) );

    // we've just left a cell -- shift index of precomputed control points
    // to the next one.
    ++cellIdx;
  }
}

Types::Coordinate
SplineWarpXform::JacobianDeterminant ( const Types::Coordinate *cp ) const
{
  const CoordinateMatrix3x3 J = this->GetJacobianAtControlPoint( cp );
  
  return this->m_InverseSpacing[0] * this->m_InverseSpacing[1] * this->m_InverseSpacing[2] * 
    ( J[0][0] * (J[1][1]*J[2][2] - J[1][2]*J[2][1]) + 
      J[0][1] * (J[1][2]*J[2][0] - J[1][0]*J[2][2]) + 
      J[0][2] * (J[1][0]*J[2][1] - J[1][1]*J[2][0]) );
}

void
SplineWarpXform
::GetJacobianConstraintThread( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t )
{
  Self::JacobianConstraintThreadInfo *info = static_cast<Self::JacobianConstraintThreadInfo*>( args );

  const SplineWarpXform *me = info->thisObject;

  const int pixelsPerRow = me->VolumeDims[0];
  std::vector<double> valuesJ( pixelsPerRow );

  const int rowCount = ( me->VolumeDims[1] * me->VolumeDims[2] );
  const int rowFrom = ( rowCount / taskCnt ) * taskIdx;
  const int rowTo = ( taskIdx == (taskCnt-1) ) ? rowCount : ( rowCount/taskCnt ) * (taskIdx+1);
  int rowsToDo = rowTo - rowFrom;

  int yFrom = rowFrom % me->VolumeDims[1];
  int zFrom = rowFrom / me->VolumeDims[2];

  double constraint = 0;
  for ( int z = zFrom; (z < me->VolumeDims[2]) && rowsToDo; ++z ) 
    {
    for ( int y = yFrom; (y < me->VolumeDims[1]) && rowsToDo; yFrom = 0, ++y, --rowsToDo ) 
      {
      me->GetJacobianDeterminantRow( &(valuesJ[0]), 0, y, z, pixelsPerRow );
      for ( int x = 0; x < pixelsPerRow; ++x ) 
	{
	constraint += fabs( log ( valuesJ[x] / me->m_GlobalScaling ) );
	}
      }
    }
  
  // Divide by number of control points to normalize with respect to the
  // number of local Jacobians in the computation.
  info->Constraint = constraint;
}

Types::Coordinate
SplineWarpXform::GetJacobianConstraint () const
{
  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  const size_t numberOfThreads = threadPool.GetNumberOfThreads();
  const size_t numberOfTasks = std::min<size_t>( 4 * numberOfThreads - 3, this->m_Dims[2] );
  
  // Info blocks for parallel tasks that evaulate the constraint.
  std::vector<Self::JacobianConstraintThreadInfo> constraintTaskInfo( numberOfTasks );
  for ( size_t taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx ) 
    {
    constraintTaskInfo[taskIdx].thisObject = this;
    }
  
  threadPool.Run( Self::GetJacobianConstraintThread, constraintTaskInfo );
  
  double constraint = 0;
  for ( size_t taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx ) 
    {
    constraint += constraintTaskInfo[taskIdx].Constraint;
    }
  
  // Divide by number of control points to normalize with respect to the
  // number of local Jacobians in the computation.
  return constraint / ( VolumeDims[0] * VolumeDims[1] * VolumeDims[2] );
}

void 
SplineWarpXform::GetJacobianConstraintDerivative
( double& lower, double& upper, const int param, const DataGrid::RegionType& voi, const Types::Coordinate step ) const
{
  const int pixelsPerRow = voi.To()[0] - voi.From()[0];
  std::vector<double> valuesJ( pixelsPerRow );
  
  double ground = 0;

  for ( int k = voi.From()[2]; k < voi.To()[2]; ++k )
    for ( int j = voi.From()[1]; j < voi.To()[1]; ++j ) 
      {
      this->GetJacobianDeterminantRow( &(valuesJ[0]), voi.From()[0], j, k, pixelsPerRow );
      for ( int i = 0; i < pixelsPerRow; ++i )
	ground += fabs( log( valuesJ[i] / this->m_GlobalScaling ) );
      }
  
  upper = -ground;
  lower = -ground;
  
  const Types::Coordinate oldCoeff = this->m_Parameters[param];
  this->m_Parameters[param] += step;
  for ( int k = voi.From()[2]; k < voi.To()[2]; ++k )
    for ( int j = voi.From()[1]; j < voi.To()[1]; ++j ) 
      {
      this->GetJacobianDeterminantRow( &(valuesJ[0]), voi.From()[0], j, k, pixelsPerRow );
      for ( int i = 0; i < pixelsPerRow; ++i )
	{
	upper += fabs( log( valuesJ[i] / this->m_GlobalScaling ) );
	}
      }
  
  this->m_Parameters[param] = oldCoeff - step;
  for ( int k = voi.From()[2]; k < voi.To()[2]; ++k )
    for ( int j = voi.From()[1]; j < voi.To()[1]; ++j ) 
      {
      this->GetJacobianDeterminantRow( &(valuesJ[0]), voi.From()[0], j, k, pixelsPerRow );
      for ( int i = 0; i < pixelsPerRow; ++i )
	{
	lower += fabs( log( valuesJ[i] / this->m_GlobalScaling ) );
	}
      }
  this->m_Parameters[param] = oldCoeff;
  
  const double invVolume = 1.0 / voi.Size();
  upper *= invVolume;
  lower *= invVolume;
}

void 
SplineWarpXform::GetJacobianConstraintDerivative
( double& lower, double& upper, const int param, const Types::Coordinate step )
  const
{
  const int controlPointIdx = param / nextI;
  const unsigned short x =  ( controlPointIdx %  this->m_Dims[0] );
  const unsigned short y = ( (controlPointIdx /  this->m_Dims[0]) % this->m_Dims[1] );
  const unsigned short z = ( (controlPointIdx /  this->m_Dims[0]) / this->m_Dims[1] );
  
  const int thisDim = param % nextI;
  const Types::Coordinate* coeff = this->m_Parameters + param - thisDim;
  
  double ground = 0;

  const int iFrom = std::max( -1, 1-x );
  const int jFrom = std::max( -1, 1-y );
  const int kFrom = std::max( -1, 1-z );

  const int iTo = std::min( 1, this->m_Dims[0]-2-x );
  const int jTo = std::min( 1, this->m_Dims[1]-2-y );
  const int kTo = std::min( 1, this->m_Dims[2]-2-z );

  for ( int k = kFrom; k < kTo; ++k )
    for ( int j = jFrom; j < jTo; ++j )
      for ( int i = iFrom; i < iTo; ++i )
	ground += fabs( log( this->GetJacobianDeterminant( Self::SpaceVectorType::FromPointer( coeff + i*nextI + j*nextJ + k*nextK ) ) / this->m_GlobalScaling ) );

  upper = -ground;
  lower = -ground;

  const Types::Coordinate oldCoeff = this->m_Parameters[param];
  this->m_Parameters[param] += step;
  for ( int k = kFrom; k < kTo; ++k )
    for ( int j = jFrom; j < jTo; ++j )
      for ( int i = iFrom; i < iTo; ++i )
	upper += fabs( log( this->GetJacobianDeterminant( Self::SpaceVectorType::FromPointer( coeff + i*nextI + j*nextJ + k*nextK ) ) / this->m_GlobalScaling ) );

  this->m_Parameters[param] = oldCoeff - step;
  for ( int k = kFrom; k < kTo; ++k )
    for ( int j = jFrom; j < jTo; ++j )
      for ( int i = iFrom; i < iTo; ++i )
	lower += fabs( log( this->GetJacobianDeterminant( Self::SpaceVectorType::FromPointer( coeff + i*nextI + j*nextJ + k*nextK ) ) / this->m_GlobalScaling ) );
  this->m_Parameters[param] = oldCoeff;

  upper /= this->m_NumberOfControlPoints;
  lower /= this->m_NumberOfControlPoints;
}

void
SplineWarpXform::RelaxToUnfold()
{
  std::vector<byte> cpList( this->m_NumberOfControlPoints );
  std::vector<double> jacobiansRow( this->VolumeDims[0] );

  bool isFolded = true;
  while ( isFolded )
    {
    // check all Jacobian determinant values to see if grid is folded, and what control points are affected
    isFolded = false;
    std::fill( cpList.begin(), cpList.end(), 0 );

    for ( int k = 0; k < this->VolumeDims[2]; ++k )
      {
      for ( int j = 0; j < this->VolumeDims[1]; ++j )
	{
	this->GetJacobianDeterminantRow( &(jacobiansRow[0]), 0, j, k, this->VolumeDims[0] );
	for ( int i = 0; i < this->VolumeDims[0]; ++i )
	  {
	  if ( jacobiansRow[i] <= 0 )
	    {
	    isFolded = true;
	    cpList[ (this->m_GridOffsets[0][i] + this->m_GridOffsets[1][j] + this->m_GridOffsets[2][k])/3 ] = 1;
	    }
	  }
	}
      }

    if ( isFolded )
      {
      // Expand affected control points to their 4x4x4 neighbourhoods on all sides.
      for ( int k = 0; k < this->m_Dims[2]; ++k )
	{
	for ( int j = 0; j < this->m_Dims[1]; ++j )
	  {
	  for ( int i = 0; i < this->m_Dims[0]; ++i )
	    {
	    if ( cpList[ (i * this->nextI + j * this->nextJ + k * this->nextK)/3 ] == 1 )
	      {
	      for ( int kk = std::max( 0, k-3 ); kk < std::min( k+4, this->m_Dims[2] ); ++kk )
		{
		for ( int jj = std::max( 0, j-3 ); jj < std::min( j+4, this->m_Dims[1] ); ++jj )
		  {
		  for ( int ii = std::max( 0, i-3 ); ii < std::min( i+4, this->m_Dims[0] ); ++ii )
		    {
		    const size_t cpIdx = (ii * this->nextI + jj * this->nextJ + kk * this->nextK)/3;
		    if ( cpList[ cpIdx ] != 1 )
		      {
		      cpList[ cpIdx ] = 2;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}

      // Get pure deformation at each control point.
      std::vector<Types::Coordinate> pureDeformation( 3 * this->m_NumberOfControlPoints );

      size_t param = 0;
      for ( int k = 0; k < this->m_Dims[2]; ++k )
	{
	for ( int j = 0; j < this->m_Dims[1]; ++j )
	  {
	  for ( int i = 0; i < this->m_Dims[0]; ++i, param+=3 )
	    {
	    const Self::SpaceVectorType cp = this->m_InitialAffineXform->GetInverse()->Apply( Self::SpaceVectorType::FromPointer( this->m_Parameters+param ) ) - this->GetOriginalControlPointPosition( i, j, k );
	    for ( int dim = 0; dim < 3; ++dim )
	      pureDeformation[param+dim] = cp[dim];
	    }
	  }
	}
      
      // regularize the affected control points
      std::vector<Types::Coordinate> smoothed( this->m_NumberOfControlPoints );
      for ( int dim = 0; dim < 3; ++dim )
	{
	// get all control point positions for the "dim" component
	size_t cp = 0;
	for ( int k = 0; k < this->m_Dims[2]; ++k )
	  {
	  for ( int j = 0; j < this->m_Dims[1]; ++j )
	    {
	    for ( int i = 0; i < this->m_Dims[0]; ++i, ++cp )
	      {
	      if ( cpList[cp] )
		{
		// if this is a control point that needs regularizing, do it
		Types::Coordinate weight = 0;
		FixedVector<3,Types::Coordinate> delta;		    
		for ( int kk = std::max( 0, k-3 ); kk < std::min( k+4, this->m_Dims[2] ); ++kk )
		  {
		  delta[2] = k-kk;
		  for ( int jj = std::max( 0, j-3 ); jj < std::min( j+4, this->m_Dims[1] ); ++jj )
		    {
		    delta[1] = j-jj;
		    for ( int ii = std::max( 0, i-3 ); ii < std::min( i+4, this->m_Dims[0] ); ++ii )
		      {
		      delta[0] = i-ii;
		      
		      const Types::Coordinate w = exp( -delta.SumOfSquares() );
		      smoothed[cp] = w * pureDeformation[ dim + ii*this->nextI + jj*this->nextJ + kk*this->nextK ];
		      weight += w;
		      }
		    }
		  }
		if ( weight > 0 )
		  smoothed[cp] /= weight;
		}
	      else
		{
		// if this control point does not need regularizing, keep it as is
		smoothed[cp] = pureDeformation[dim+cp*3];
		}
	      }
	    }
	  }
	
	// copy modified control point position component back
	for ( size_t cp = 0; cp < this->m_NumberOfControlPoints; ++cp )
	  {
	  pureDeformation[dim+cp*3] = smoothed[cp];
	  }
	}
      
      // put offsets and affine component back
      param = 0;
      for ( int k = 0; k < this->m_Dims[2]; ++k )
	{
	for ( int j = 0; j < this->m_Dims[1]; ++j )
	  {
	  for ( int i = 0; i < this->m_Dims[0]; ++i, param+=3 )
	    {
	    Self::SpaceVectorType cp = Self::SpaceVectorType::FromPointer( &(pureDeformation[0])+param );
	    cp += this->GetOriginalControlPointPosition( i, j, k );
	    cp = this->m_InitialAffineXform->Apply( cp );
	    
	    for ( int dim = 0; dim < 3; ++dim )
	      this->m_Parameters[param+dim] = cp[dim];
	    }
	  }
	}
      }
    }
}

} // namespace cmtk
