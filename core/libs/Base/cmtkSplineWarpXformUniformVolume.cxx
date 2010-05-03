/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include <cmtkSplineWarpXformUniformVolume.h>

cmtk::SplineWarpXformUniformVolume::SplineWarpXformUniformVolume( const UniformVolume& volume, const SplineWarpXform::SmartPtr& xform )
  : m_Xform( xform )
{
  this->RegisterVolume( volume );
}

void
cmtk::SplineWarpXformUniformVolume::RegisterVolumeAxis 
( const int dim, const Types::Coordinate delta, const Types::Coordinate origin, const int cpgDim, const Types::Coordinate invCpgSpacing,
  std::vector<int>& g, std::vector<Types::Coordinate>& spline, std::vector<Types::Coordinate>& dspline )
{
  g.resize( dim+1 );
  spline.resize( 4*dim );
  dspline.resize( 4*dim );

  for ( int idx=0; idx<dim; ++idx ) 
    {
    Types::Coordinate r = invCpgSpacing * (origin + delta * idx);
    g[idx] = std::min( static_cast<int>( r ), cpgDim-4 );
    const Types::Coordinate f = r - g[idx];
    for ( int k = 0; k < 4; ++k ) 
      {
      spline[4*idx+k] = CubicSpline::ApproxSpline( k, f );
      dspline[4*idx+k] = CubicSpline::DerivApproxSpline( k, f );
      }
    }
  // guard element
  g[dim] = -1;
}

void
cmtk::SplineWarpXformUniformVolume::RegisterVolume
( const UniformVolume& volume )
{
  const SplineWarpXform& xform = *(this->m_Xform);

  this->RegisterVolumeAxis( volume.m_Dims[0], volume.m_Delta[0], volume.m_Offset[0], xform.m_Dims[0], xform.InverseSpacing[0], this->gX, this->splineX, this->dsplineX );
  this->RegisterVolumeAxis( volume.m_Dims[1], volume.m_Delta[1], volume.m_Offset[1], xform.m_Dims[1], xform.InverseSpacing[1], this->gY, this->splineY, this->dsplineY );
  this->RegisterVolumeAxis( volume.m_Dims[2], volume.m_Delta[2], volume.m_Offset[2], xform.m_Dims[2], xform.InverseSpacing[2], this->gZ, this->splineZ, this->dsplineZ );
  
  for ( int idx = 0; idx < volume.m_Dims[0]; ++idx ) 
    gX[idx] *= xform.nextI;
  for ( int idx = 0; idx < volume.m_Dims[1]; ++idx ) 
    gY[idx] *= xform.nextJ;
  for ( int idx = 0; idx < volume.m_Dims[2]; ++idx ) 
    gZ[idx] *= xform.nextK;
}

void
cmtk::SplineWarpXformUniformVolume
::GetTransformedGrid 
( Vector3D& v, const int idxX, const int idxY, const int idxZ ) const
{
  const SplineWarpXform& xform = *(this->m_Xform);

  const Types::Coordinate* coeff = xform.m_Parameters + gX[idxX] + gY[idxY] + gZ[idxZ];
  const Types::Coordinate *spX = &splineX[idxX<<2], *spY = &splineY[idxY<<2], *spZ = &splineZ[idxZ<<2];
  
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
	  coeff_ll += xform.nextJ;
	  }	
	mm += spZ[m] * ll;
	coeff_mm += xform.nextK;
	}
      v.XYZ[ dim ] = mm;
      ++coeff;
    }
}

void 
cmtk::SplineWarpXformUniformVolume::GetTransformedGridSequence
( Vector3D *const vIn, const size_t numPoints, const int idxX, const int idxY, const int idxZ ) 
  const
{
  Vector3D *v = vIn;

  const SplineWarpXform& xform = *(this->m_Xform);
  const Types::Coordinate* coeff = xform.m_Parameters + gX[idxX] + gY[idxY] + gZ[idxZ];
  const Types::Coordinate *spX = &splineX[idxX<<2], *spY = &splineY[idxY<<2], *spZ = &splineZ[idxZ<<2];
  
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
  const int numberOfCells = (gX[idxX + numPoints - 1] - gX[idxX]) / xform.nextI + 4;
  
  // pre-compute the contributions of all control points in y- and z-direction
  // along the way
  Types::Coordinate phiComp;
  std::vector<Types::Coordinate> phiHat( 3*numberOfCells );

  const int *gpo;
  int phiIdx = 0;
  for ( int cell = 0; cell < numberOfCells; ++cell, coeff += xform.nextI ) 
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
      Types::Coordinate* vPtr = v->XYZ;
      // compute transformed voxel by taking precomputed y- and z-contributions
      // and adding x. The loops to do this have been unrolled for increased
      // performance.
      vPtr[0] = spX[0] * phiPtr[0] + spX[1] * phiPtr[3] + spX[2] * phiPtr[6] + spX[3] * phiPtr[9];
      vPtr[1] = spX[0] * phiPtr[1] + spX[1] * phiPtr[4] + spX[2] * phiPtr[7] + spX[3] * phiPtr[10];
      vPtr[2] = spX[0] * phiPtr[2] + spX[1] * phiPtr[5] + spX[2] * phiPtr[8] + spX[3] * phiPtr[11];
      
      // go to next voxel
      ++i;
      spX += 4;
      ++v;
      // repeat this until we leave current CPG cell.
      } 
    while ( ( gX[i-1] == gX[i] ) && ( i < lastPoint ) );
    
    // we've just left a cell -- shift index of precomputed control points
    // to the next one.
    ++cellIdx;
    }
}

