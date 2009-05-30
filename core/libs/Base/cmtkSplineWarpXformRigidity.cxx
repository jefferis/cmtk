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

#include <cmtkArray.h>
#include <cmtkQRDecomposition.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

Types::Coordinate
SplineWarpXform::GetRigidityConstraintSparse () const
{
  double Constraint = 0;
  CoordinateMatrix3x3 J;

  const Types::Coordinate* coeff = this->m_Parameters + nextI + nextJ + nextK;
  for ( int z = 1; z<Dims[2]-1; ++z, coeff+=2*nextJ )
    for ( int y = 1; y<Dims[1]-1; ++y, coeff+=2*nextI )
      for ( int x = 1; x<Dims[0]-1; ++x, coeff+=nextI )
	{
	this->GetJacobian( coeff, J );
	Constraint += this->GetRigidityConstraint( J );
	}
  
  // Divide by number of control points to normalize with respect to the
  // number of local Jacobians in the computation.
  return (Types::Coordinate)(Constraint / this->NumberOfControlPoints);
}

void 
SplineWarpXform::GetRigidityConstraintDerivative
( double& lower, double& upper, const int param, const Rect3D& voi,
  const Types::Coordinate step ) const
{
  const int pixelsPerRow = voi.endX - voi.startX;
  Array<CoordinateMatrix3x3> arrayJ( pixelsPerRow );
  
  double ground = 0;

  for ( int k = voi.startZ; k < voi.endZ; ++k )
    for ( int j = voi.startY; j < voi.endY; ++j ) 
      {
      this->GetJacobianSequence( arrayJ, voi.startX, j, k, pixelsPerRow );
      for ( int i = 0; i < pixelsPerRow; ++i ) 
	{
	ground += this->GetRigidityConstraint( arrayJ[i] );
	}
      }
  
  upper = -ground;
  lower = -ground;
  
  const Types::Coordinate oldCoeff = this->m_Parameters[param];
  this->m_Parameters[param] += step;
  for ( int k = voi.startZ; k < voi.endZ; ++k )
    for ( int j = voi.startY; j < voi.endY; ++j ) 
      {
      this->GetJacobianSequence( arrayJ, voi.startX, j, k, pixelsPerRow );
      for ( int i = 0; i < pixelsPerRow; ++i ) 
	{
	upper += this->GetRigidityConstraint( arrayJ[i] );
	}
      }
  
  this->m_Parameters[param] = oldCoeff - step;
  for ( int k = voi.startZ; k < voi.endZ; ++k )
    for ( int j = voi.startY; j < voi.endY; ++j ) 
      {
      this->GetJacobianSequence( arrayJ, voi.startX, j, k, pixelsPerRow );
      for ( int i = 0; i < pixelsPerRow; ++i ) 
	{
	lower += this->GetRigidityConstraint( arrayJ[i] );
	}
      }
  this->m_Parameters[param] = oldCoeff;

  const double invVolume = 1.0 / ((voi.endX-voi.startX)*(voi.endY-voi.startY)*(voi.endZ-voi.startZ));
  upper *= invVolume;
  lower *= invVolume;
}

void 
SplineWarpXform::GetRigidityConstraintDerivative
( double& lower, double& upper, const int param, const Rect3D& voi,
  const Types::Coordinate step, const DataGrid* weightMap ) const
{
  const int pixelsPerRow = voi.endX - voi.startX;
  Array<CoordinateMatrix3x3> arrayJ( pixelsPerRow );
  
  double ground = 0;

  for ( int k = voi.startZ; k < voi.endZ; ++k )
    for ( int j = voi.startY; j < voi.endY; ++j ) 
      {
      this->GetJacobianSequence( arrayJ, voi.startX, j, k, pixelsPerRow );
      for ( int i = 0; i < pixelsPerRow; ++i ) 
	{
	ground += weightMap->GetDataAt( voi.startX + i, j, k ) * this->GetRigidityConstraint( arrayJ[i] );
	}
      }
  
  upper = -ground;
  lower = -ground;
  
  const Types::Coordinate oldCoeff = this->m_Parameters[param];
  this->m_Parameters[param] += step;
  for ( int k = voi.startZ; k < voi.endZ; ++k )
    for ( int j = voi.startY; j < voi.endY; ++j ) 
      {
      this->GetJacobianSequence( arrayJ, voi.startX, j, k, pixelsPerRow );
      for ( int i = 0; i < pixelsPerRow; ++i ) 
	{
	upper += weightMap->GetDataAt( voi.startX + i, j, k ) *
	  this->GetRigidityConstraint( arrayJ[i] );
	}
      }
  
  this->m_Parameters[param] = oldCoeff - step;
  for ( int k = voi.startZ; k < voi.endZ; ++k )
    for ( int j = voi.startY; j < voi.endY; ++j ) 
      {
      this->GetJacobianSequence( arrayJ, voi.startX, j, k, pixelsPerRow );
      for ( int i = 0; i < pixelsPerRow; ++i ) 
	{
	lower += weightMap->GetDataAt( voi.startX + i, j, k ) * this->GetRigidityConstraint( arrayJ[i] );
	}
      }
  this->m_Parameters[param] = oldCoeff;

  double invVolume = 1.0 / ((voi.endX-voi.startX)*(voi.endY-voi.startY)*(voi.endZ-voi.startZ));
  upper *= invVolume;
  lower *= invVolume;
}

void 
SplineWarpXform::GetRigidityConstraintDerivative
( double& lower, double& upper, const int param, const Types::Coordinate step )
  const
{
  const int controlPointIdx = param / nextI;
  const unsigned short x =  ( controlPointIdx %  Dims[0] );
  const unsigned short y = ( (controlPointIdx /  Dims[0]) % Dims[1] );
  const unsigned short z = ( (controlPointIdx /  Dims[0]) / Dims[1] );
  
  const int thisDim = param % nextI;
  const Types::Coordinate* coeff = this->m_Parameters + param - thisDim;
  
  double ground = 0;

  const int iFrom = std::max( -1, 1-x );
  const int jFrom = std::max( -1, 1-y );
  const int kFrom = std::max( -1, 1-z );

  const int iTo = std::min( 1, Dims[0]-2-x );
  const int jTo = std::min( 1, Dims[1]-2-y );
  const int kTo = std::min( 1, Dims[2]-2-z );

  CoordinateMatrix3x3 J;
  for ( int k = kFrom; k < kTo; ++k )
    for ( int j = jFrom; j < jTo; ++j )
      for ( int i = iFrom; i < iTo; ++i )
	{
	this->GetJacobianAtControlPoint( coeff + i*nextI + j*nextJ + k*nextK, J );
	ground += this->GetRigidityConstraint( J );
	}
  
  upper = -ground;
  lower = -ground;
  
  const Types::Coordinate oldCoeff = this->m_Parameters[param];
  this->m_Parameters[param] += step;
  for ( int k = kFrom; k < kTo; ++k )
    for ( int j = jFrom; j < jTo; ++j )
      for ( int i = iFrom; i < iTo; ++i )
	{
	this->GetJacobianAtControlPoint( coeff + i*nextI + j*nextJ + k*nextK, J );
	upper += this->GetRigidityConstraint( J );
	}

  this->m_Parameters[param] = oldCoeff - step;
  for ( int k = kFrom; k < kTo; ++k )
    for ( int j = jFrom; j < jTo; ++j )
      for ( int i = iFrom; i < iTo; ++i )
	{
	this->GetJacobianAtControlPoint( coeff + i*nextI + j*nextJ + k*nextK, J );
	lower += this->GetRigidityConstraint( J );
	}
  this->m_Parameters[param] = oldCoeff;

  upper /= this->NumberOfControlPoints;
  lower /= this->NumberOfControlPoints;
}

Types::Coordinate
SplineWarpXform::GetRigidityConstraint () const
{
  const int pixelsPerRow = VolumeDims[0];
  Array<CoordinateMatrix3x3> arrayJ( pixelsPerRow );

  double constraint = 0;
  for ( int z = 0; z < VolumeDims[2]; ++z )
    for ( int y = 0; y < VolumeDims[1]; ++y ) 
      {
      this->GetJacobianSequence( arrayJ, 0, y, z, pixelsPerRow );
      for ( int x = 0; x < pixelsPerRow; ++x ) 
	{
	constraint += this->GetRigidityConstraint( arrayJ[x] );
	}
      }
  
  // Divide by number of control points to normalize with respect to the
  // number of local Jacobians in the computation.
  return constraint / ( VolumeDims[0] * VolumeDims[1] * VolumeDims[2] );
}

Types::Coordinate
SplineWarpXform::GetRigidityConstraint( const DataGrid* weightMap ) const
{
  const int pixelsPerRow = VolumeDims[0];
  Array<CoordinateMatrix3x3> arrayJ( pixelsPerRow );

  double constraint = 0;
  for ( int z = 0; z < VolumeDims[2]; ++z )
    for ( int y = 0; y < VolumeDims[1]; ++y ) 
      {
      this->GetJacobianSequence( arrayJ, 0, y, z, pixelsPerRow );
      for ( int x = 0; x < pixelsPerRow; ++x ) 
	{
	constraint += weightMap->GetDataAt( x, y, z ) * this->GetRigidityConstraint( arrayJ[x] );
	}
      }
  
  // Divide by number of control points to normalize with respect to the
  // number of local Jacobians in the computation.
  return constraint / ( VolumeDims[0] * VolumeDims[1] * VolumeDims[2] );
}

Types::Coordinate
SplineWarpXform::GetRigidityConstraint( const CoordinateMatrix3x3& J ) 
  const
{
  Matrix2D<Types::Coordinate> matrix2d( 3, 3 );
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      matrix2d[i][j] = J[i][j];

  QRDecomposition<Types::Coordinate> qr( matrix2d );
  Matrix2D<Types::Coordinate> R = *(qr.GetR());
  
  Array<Types::Coordinate> R_diagonal( 3 );
  for ( int i = 0; i < 3; i++ )
    R_diagonal[i] = R[i][i];

  return MathUtil::Square( R[0][1] / R_diagonal[0] ) + MathUtil::Square( R[0][2] / R_diagonal[0] ) + MathUtil::Square( R[1][2] / R_diagonal[1] );
}

} // namespace cmtk
