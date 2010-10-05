/*
//
//  Copyright 2010 SRI International
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

#include <Base/cmtkParametricPlane.h>

#include <algorithm>
#include <math.h>
#include <iostream>

// test parametric plane reflection with zero offset
int
testParametricPlaneMirror()
{
  typedef cmtk::ParametricPlane PlaneType;
  typedef PlaneType::CoordinateVectorType VectorType;
  typedef cmtk::AffineXform::MatrixType MatrixType;

  VectorType origin;
  std::fill( origin.begin(), origin.end(), 0 );

  PlaneType plane;
  plane.SetOrigin( origin );
  plane.SetTheta( cmtk::Units::Degrees( 30 ) );
  plane.SetPhi( cmtk::Units::Degrees( 45 ) );

  for ( double rho = 0; rho < 21; rho += 20 )
    {
    plane.SetRho( rho );

    VectorType v;
    v[0] = v[1] = v[2] = 1;
    
    VectorType pv = v;
    plane.MirrorInPlace( pv );
    
    MatrixType pm = plane.GetMirrorXformMatrix();
    
    VectorType mv = v * pm;
    
    if ( sqrt( (mv-pv).SumOfSquares() ) > 1e-5 )
      {
      std::cerr << "ParametricPlane::MirrorInPlace and mirror matrix produce different results." << std::endl;
      return 1;
      }
    
    mv *= pm;
    if ( sqrt( (mv-v).SumOfSquares() ) > 1e-5 )
      {
      std::cerr << "ParametricPlane mirror matrix applied twice does not return to original point." << std::endl;
      return 1;
      }
    
    plane.MirrorInPlace( pv );
    if ( sqrt( (pv-v).SumOfSquares() ) > 1e-5 )
      {
      std::cerr << "ParametricPlane::MirrorInPlane applied twice does not return to original point." << std::endl;
      return 1;
      }
    }
  
  return 0;
}

// test parametric plane reflection with non-zero offset
int
testParametricPlaneMirrorOffset()
{
  typedef cmtk::ParametricPlane PlaneType;
  typedef PlaneType::CoordinateVectorType VectorType;
  typedef cmtk::AffineXform::MatrixType MatrixType;

  VectorType origin;
  origin[0] = 10;
  origin[1] = 20;
  origin[2] = 30;

  PlaneType plane;
  plane.SetOrigin( origin );
  plane.SetTheta( cmtk::Units::Degrees( 30 ) );
  plane.SetPhi( cmtk::Units::Degrees( 45 ) );

  for ( double rho = 0; rho < 21; rho += 20 )
    {
    plane.SetRho( rho );

    VectorType v;
    v[0] = v[1] = v[2] = 1;
    
    VectorType pv = v;
    plane.MirrorInPlace( pv );
    
    MatrixType pm = plane.GetMirrorXformMatrix();
    
    VectorType mv = v * pm;
    
    if ( sqrt( (mv-pv).SumOfSquares() ) > 1e-5 )
      {
      std::cerr << "ParametricPlane::MirrorInPlace and mirror matrix produce different results." << std::endl;
      return 1;
      }
    
    mv *= pm;
    if ( sqrt( (mv-v).SumOfSquares() ) > 1e-5 )
      {
      std::cerr << "ParametricPlane mirror matrix applied twice does not return to original point." << std::endl;
      return 1;
      }
    
    plane.MirrorInPlace( pv );
    if ( sqrt( (pv-v).SumOfSquares() ) > 1e-5 )
      {
      std::cerr << "ParametricPlane::MirrorInPlane applied twice does not return to original point." << std::endl;
      return 1;
      }
    }
  
  return 0;
}

