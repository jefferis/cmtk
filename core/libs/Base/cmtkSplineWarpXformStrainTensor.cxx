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

#include <cmtkSplineWarpXform.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

CoordinateMatrix3x3 
SplineWarpXform::GetStrainTensor( const Vector3D& v ) const
{
  CoordinateMatrix3x3 J;
  this->GetStrainTensor( v, J );
  return J;
}

void
SplineWarpXform::GetStrainTensor( const Vector3D& v, CoordinateMatrix3x3& t ) const
{
  Types::Coordinate r[3], f[3];
  int grid[3];
  
  for ( int dim = 0; dim<3; ++dim ) 
    {
    r[dim] = this->InverseSpacing[dim] * v.XYZ[dim];
    grid[dim] = std::min( static_cast<int>( r[dim] ), Dims[dim]-4 );
    f[dim] = r[dim] - grid[dim];
    assert( (f[dim] >= 0.0) && (f[dim] <= 1.0) );
    }

  const Types::Coordinate* coeff = 
    this->m_Parameters + 3 * ( grid[0] + Dims[0] * (grid[1] + Dims[1] * grid[2]) );

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
	  kk[1] += CubicSpline::ApproxSpline( k, f[0] ) * (*coeff_kk);
	  // dT / dz
	  kk[2] += CubicSpline::ApproxSpline( k, f[0] ) * (*coeff_kk);
	  }
	
	// dT / dx
	ll[0] += CubicSpline::ApproxSpline( l, f[1] ) * kk[0];
	// dT / dy
	ll[1] += CubicSpline::DerivApproxSpline( l, f[1] ) * kk[1];
	// dT / dz
	ll[2] += CubicSpline::ApproxSpline( l, f[1] ) * kk[2];
	}

      // dT / dx
      t[dim][0] += CubicSpline::ApproxSpline( m, f[2] ) * ll[0];
      // dT / dy
      t[dim][1] += CubicSpline::ApproxSpline( m, f[2] ) * ll[1];
      // dT / dz
      t[dim][2] += CubicSpline::DerivApproxSpline( m, f[2] ) * ll[2];
      }
    }

  // scale with grid spacing to normalize strain tensor (chain rule of derivation)
  // and compute J-I
  for ( int i = 0; i<3; ++i ) 
    {
    for ( int j = 0; j<3; ++j )
      t[i][j] *= this->InverseSpacing[j];

    // subtract diagonal to eliminate constant component
    t[i][i] -= 1.0;
    }

  // compute 0.5 * (J-I)+(J-I)^T
  for ( int i = 0; i<3; ++i ) 
    {
    for ( int j = 0; j<i; ++j )
      {
      t[i][j] += t[j][i];
      t[i][j] *= 0.5;
      }
    t[i][i] *= 0.5;
    }

  for ( int i = 0; i<3; ++i ) 
    {
    for ( int j = 0; j<i; ++j )
      t[j][i] = t[i][j];
    }
}

} // namespace cmtk
