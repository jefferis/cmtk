/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include "cmtkDeformationField.h"

#include <Base/cmtkCubicSpline.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
DeformationField::GetJacobian
( const Self::SpaceVectorType& v, CoordinateMatrix3x3& J ) const
{
  Types::Coordinate r[3], f[3];
  int grid[3];
  
  for ( int dim = 0; dim<3; ++dim ) 
    {
    r[dim] = this->m_InverseSpacing[dim] * (v[dim] - this->m_Offset[dim]);
    grid[dim] = static_cast<int>( r[dim]-1 );
    if ( (grid[dim] < 0) || (grid[dim] >= this->m_Dims[dim]-3) )
      {
      J.Fill( 0.0 );
      return;
      }
    f[dim] = r[dim] - grid[dim] - 1;
    }

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
	  kk[0] += CubicSpline::DerivInterpSpline( k, f[0] ) * (*coeff_kk);
	  // dT / dy
	  const Types::Coordinate tmp = CubicSpline::InterpSpline( k, f[0] ) * (*coeff_kk);
	  kk[1] += tmp;
	  // dT / dz
	  kk[2] += tmp;
	  }
	
	// dT / dx
	const Types::Coordinate tmp = CubicSpline::InterpSpline( l, f[1] );
	ll[0] += tmp * kk[0];
	// dT / dy
	ll[1] += CubicSpline::DerivInterpSpline( l, f[1] ) * kk[1];
	// dT / dz
	ll[2] += tmp * kk[2];
	}

      // dT / dx
      const Types::Coordinate tmp = CubicSpline::InterpSpline( m, f[2] );
      J[dim][0] += tmp * ll[0];
      // dT / dy
      J[dim][1] += tmp * ll[1];
      // dT / dz
      J[dim][2] += CubicSpline::DerivInterpSpline( m, f[2] ) * ll[2];
      }
    }
}

} // namespace cmtk
