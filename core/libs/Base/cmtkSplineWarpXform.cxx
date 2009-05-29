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

#include <cmtkSplineWarpXform.h>

#include <cmtkMathUtil.h>
#include <cmtkBitVector.h>
#include <cmtkArray.h>

#include <assert.h>

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
  this->GlobalScaling = 1.0;
}

SplineWarpXform::SplineWarpXform
( const UniformVolume *volume, const Types::Coordinate delta, const bool exactDelta )
{
  this->Init();
  memcpy( Domain, volume->Size, sizeof( Domain ) );
  InitialAffineXform = AffineXform::SmartPtr( NULL );

  if ( exactDelta ) 
    {
    for ( int dim=0; dim<3; ++dim ) 
      {
      Spacing[dim] = delta;
      Dims[dim] = static_cast<int>( 4 + (Domain[dim] / Spacing[dim]) );
      Domain[dim] = (Dims[dim] - 3) * Spacing[dim];
      }
    } 
  else
    {
    for ( int dim=0; dim<3; ++dim )
      Dims[dim] = 2 + std::max( 2, 1+static_cast<int>( Domain[dim]/delta ) );
    }
  
  NumberOfControlPoints = Dims[0] * Dims[1] * Dims[2];
  this->AllocateParameterVector( 3 * NumberOfControlPoints );
  
  this->Update( exactDelta );
  this->InitControlPoints( InitialAffineXform );
  this->RegisterVolume( volume );
}

SplineWarpXform::SplineWarpXform 
( const Types::Coordinate domain[3], const Types::Coordinate delta,
  const AffineXform* initialXform, const bool exactDelta  )
{
  this->Init( domain, delta, initialXform, exactDelta );
}

void
SplineWarpXform
::Init
( const Types::Coordinate domain[3], const Types::Coordinate delta, const AffineXform* initialXform, const bool exactDelta  )
{
  this->Init();
  memcpy( Domain, domain, sizeof( Domain ) );
  InitialAffineXform = AffineXform::SmartPtr( initialXform->Clone() );

  if ( exactDelta ) 
    {
    for ( int dim=0; dim<3; ++dim ) 
      {
      Spacing[dim] = delta;
      Dims[dim] = static_cast<int>( 4 + (Domain[dim] / Spacing[dim]) );
      Domain[dim] = (Dims[dim] - 3) * Spacing[dim];
      }
    } 
  else
    {
    for ( int dim=0; dim<3; ++dim )
      Dims[dim] = 2 + std::max( 2, 1+static_cast<int>( domain[dim]/delta ) );
    }
  
  NumberOfControlPoints = Dims[0] * Dims[1] * Dims[2];
  this->AllocateParameterVector( 3 * NumberOfControlPoints );
  
  this->Update( exactDelta );
  this->InitControlPoints( InitialAffineXform );
}

SplineWarpXform::SplineWarpXform
( const Types::Coordinate domain[3], const int dims[3],
  CoordinateVector::SmartPtr& parameters, 
  const AffineXform* initialXform )
{
  this->Init();
  memcpy( Domain, domain, sizeof( Domain ) );
  memcpy( Dims, dims, sizeof( Dims ) );

  if ( initialXform )
    {
    InitialAffineXform = AffineXform::SmartPtr( initialXform->Clone() );
    GlobalScaling = InitialAffineXform->GetGlobalScaling();
    } 
  else
    {
    InitialAffineXform = AffineXform::SmartPtr( NULL );
    }

  NumberOfControlPoints = Dims[0] * Dims[1] * Dims[2];
  this->m_NumberOfParameters = 3 * NumberOfControlPoints;

  if ( parameters.IsNull() )
    this->m_ParameterVector = CoordinateVector::SmartPtr( new CoordinateVector( this->m_NumberOfParameters ) );
  else
    this->m_ParameterVector = parameters;
  this->m_Parameters = this->m_ParameterVector->Elements;

  this->Update( false /* exactDelta */ );

  if ( parameters.IsNull() )
      this->InitControlPoints( InitialAffineXform );
}

void
SplineWarpXform::InitControlPoints( const AffineXform* affineXform )
{
  Types::Coordinate *ofs = this->m_Parameters;
  Types::Coordinate pZ = -Spacing[2];
  for ( int z=0; z<Dims[2]; ++z, pZ+=Spacing[2] ) 
    {
    Types::Coordinate pY = -Spacing[1];
    for ( int y=0; y<Dims[1]; ++y, pY+=Spacing[1] ) 
      {
      Types::Coordinate pX = -Spacing[0];
      for ( int x=0; x<Dims[0]; ++x, pX+=Spacing[0], ofs+=3 ) 
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
    for ( unsigned int idx = 0; idx < NumberOfControlPoints; ++idx, ofs+=3 ) 
      {
      Vector3D p( ofs );
      affineXform->ApplyInPlaceNonVirtual( p );
      ofs[0] = p[0];
      ofs[1] = p[1];
      ofs[2] = p[2];
      }
    
    affineXform->GetScales( this->InverseAffineScaling );
    GlobalScaling = affineXform->GetGlobalScaling();
    } 
  else
    {
    InverseAffineScaling[0] = InverseAffineScaling[1] = InverseAffineScaling[2] = GlobalScaling = 1.0;
    }
}

void
SplineWarpXform::Update 
( const bool exactDelta ) 
{
  this->WarpXform::Update();

  for ( int dim=0; dim<3; ++dim ) 
    {
    assert( Dims[dim] > 3 );
    if ( exactDelta ) 
      {
      InverseSpacing[dim] = 1.0 / Spacing[dim];
      } 
    else
      {
      Spacing[dim] = Domain[dim] / (Dims[dim]-3);
      InverseSpacing[dim] = 1.0*(Dims[dim]-3) / Domain[dim];
      }
    m_Origin[dim] = -Spacing[dim];
    }
  
  int dml = 0;
  for ( int dim = 0; dim<3; ++dim )
    for ( int m = 0; m < 4; ++m )
      for ( int l = 0; l < 4; ++l, ++dml )
	GridPointOffset[dml] = dim + l * nextJ + m * nextK;
}

SplineWarpXform* 
SplineWarpXform::Clone () const
{
  SplineWarpXform *newXform = new SplineWarpXform();

  newXform->m_ParameterVector = CoordinateVector::SmartPtr( this->m_ParameterVector->Clone() );
  newXform->m_Parameters = newXform->m_ParameterVector->Elements;
  newXform->m_NumberOfParameters = this->m_NumberOfParameters;
  newXform->NumberOfControlPoints = this->NumberOfControlPoints;
  
  memcpy( newXform->Dims, Dims, sizeof( newXform->Dims ) );
  memcpy( newXform->Domain, Domain, sizeof( newXform->Domain ) );
  memcpy( newXform->Spacing, Spacing, sizeof( newXform->Spacing ) );
  memcpy( newXform->InverseSpacing, InverseSpacing, sizeof( newXform->InverseSpacing ) );
  newXform->m_Origin = this->m_Origin;

  if ( ActiveFlags ) 
    {
    BitVector::SmartPtr activeFlags( this->ActiveFlags->Clone() );
    newXform->SetActiveFlags( activeFlags );
    }
  newXform->IgnoreEdge = this->IgnoreEdge;
  newXform->FastMode = this->FastMode;

  if ( this->InitialAffineXform ) 
    {
    newXform->InitialAffineXform = AffineXform::SmartPtr( this->InitialAffineXform->Clone() );
    }
  
  newXform->GlobalScaling = this->GlobalScaling;

  newXform->nextI = this->nextI;
  newXform->nextJ = this->nextJ;
  newXform->nextK = this->nextK;
  newXform->nextIJ = this->nextIJ;
  newXform->nextIK = this->nextIK;
  newXform->nextJK = this->nextJK;
  newXform->nextIJK = this->nextIJK;
  memcpy( newXform->GridPointOffset, this->GridPointOffset, sizeof( this->GridPointOffset ) );

  memcpy( newXform->VolumeDims, this->VolumeDims, sizeof( this->VolumeDims ) );

  newXform->gX = this->gX;
  newXform->gY = this->gY;
  newXform->gZ = this->gZ;

  newXform->splineX = this->splineX;
  newXform->splineY = this->splineY;
  newXform->splineZ = this->splineZ;

  newXform->dsplineX = this->dsplineX;
  newXform->dsplineY = this->dsplineY;
  newXform->dsplineZ = this->dsplineZ;

  return newXform;
}

void
SplineWarpXform::Refine ( const int factor )
{
  if ( !this->m_ParameterVector ) return;

  if ( factor != 2 )
    {
    fputs( "WARNING: Cannot refine spline warps by factors other than 2.\n",
	   stderr );
    return;
    }
  
  int newDims[3];
  for ( int dim=0; dim<3; ++dim ) 
    newDims[dim] = 2 * Dims[dim] - 3;

  const int newNumSamples = newDims[0] * newDims[1] * newDims[2];
  const int newNumCoefficients = 3 * newNumSamples;

  CoordinateVector::SmartPtr newParameters( new CoordinateVector( newNumCoefficients ) );
  Types::Coordinate* newCoefficients = newParameters->Elements;

  Types::Coordinate newSpacing[3];
  for ( int dim=0; dim<3; ++dim ) 
    {
    newSpacing[dim] = Domain[dim] / (newDims[dim]-3);
    }

  // no linear refinement here
  const int newNextI = 3;
  const int newNextJ = newNextI * newDims[0];
  const int newNextK = newNextJ * newDims[1];
  const int newNextIJ = newNextJ + newNextI;
  const int newNextIK = newNextK + newNextI;
  const int newNextJK = newNextK + newNextJ;
  const int newNextIJK = newNextJK + newNextI;

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
	  Types::Coordinate level0[3][3];
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
	  
	  Types::Coordinate level1[3];
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

  NumberOfControlPoints = newNumSamples;
  this->m_NumberOfParameters = newNumCoefficients;

  this->m_ParameterVector = newParameters;
  this->m_Parameters = newCoefficients;

  this->DeleteParameterActiveFlags();
  memcpy( Dims, newDims, sizeof(Dims) );

  for ( int dim=0; dim<3; ++dim ) 
    {
    assert( Dims[dim] > 1 );
    Spacing[dim] = newSpacing[dim];
    InverseSpacing[dim] = 1.0 / Spacing[dim];
    m_Origin[dim] = -Spacing[dim];
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

  if ( IgnoreEdge )
    IgnoreEdge = 2 * IgnoreEdge - 1;

  this->UnRegisterVolume();
}

void 
SplineWarpXform::GetVolumeOfInfluence 
( const size_t idx, const Vector3D& fromVol, const Vector3D& toVol,
  Vector3D& fromVOI, Vector3D& toVOI, const int fastMode ) const
{
  int relIdx = idx / 3;

  int x =  ( relIdx %  Dims[0] );
  int y = ( (relIdx /  Dims[0]) % Dims[1] );
  int z = ( (relIdx /  Dims[0]) / Dims[1] );

  Types::Coordinate xLow, yLow, zLow, xUp, yUp, zUp;

  if ( (fastMode==1) || (FastMode && (fastMode<0)) ) 
    {
    xLow = Spacing[0] * std::max( 0, x-2 );
    yLow = Spacing[1] * std::max( 0, y-2 );
    zLow = Spacing[2] * std::max( 0, z-2 );
    
    xUp = Spacing[0] * std::min( Dims[0]-3, x );
    yUp = Spacing[1] * std::min( Dims[1]-3, y );
    zUp = Spacing[2] * std::min( Dims[2]-3, z );
    } 
  else
    {
    xLow = Spacing[0] * std::max( 0, x-3 );
    yLow = Spacing[1] * std::max( 0, y-3 );
    zLow = Spacing[2] * std::max( 0, z-3 );
    
    xUp = Spacing[0] * std::min( Dims[0]-2, x+1 );
    yUp = Spacing[1] * std::min( Dims[1]-2, y+1 );
    zUp = Spacing[2] * std::min( Dims[2]-2, z+1 );
    }
  
  fromVOI.Set
    ( std::min( toVol[0], std::max(xLow, fromVol[0]) ), 
      std::min( toVol[1], std::max(yLow, fromVol[1]) ), 
      std::min( toVol[2], std::max(zLow, fromVol[2]) ) );
  toVOI.Set
    ( std::max( fromVol[0], std::min(xUp, toVol[0]) ), 
      std::max( fromVol[1], std::min(yUp, toVol[1]) ), 
      std::max( fromVol[2], std::min(zUp, toVol[2]) ) );
}

void
SplineWarpXform::RegisterVolumeAxis 
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
SplineWarpXform::RegisterVolumePoints
( const int volDims[3], const Types::Coordinate delta[3], const Types::Coordinate origin[3] )
{
  this->RegisterVolumeAxis( volDims[0], delta[0], origin[0], Dims[0], InverseSpacing[0], gX, splineX, dsplineX );
  this->RegisterVolumeAxis( volDims[1], delta[1], origin[1], Dims[1], InverseSpacing[1], gY, splineY, dsplineY );
  this->RegisterVolumeAxis( volDims[2], delta[2], origin[2], Dims[2], InverseSpacing[2], gZ, splineZ, dsplineZ );

  for ( int idx = 0; idx<volDims[0]; ++idx ) gX[idx] *= nextI;
  for ( int idx = 0; idx<volDims[1]; ++idx ) gY[idx] *= nextJ;
  for ( int idx = 0; idx<volDims[2]; ++idx ) gZ[idx] *= nextK;

  memcpy( VolumeDims, volDims, sizeof( VolumeDims ) );
}

void SplineWarpXform::UnRegisterVolume()
{
  gX.resize( 0 );
  gY.resize( 0 );
  gZ.resize( 0 );

  splineX.resize( 0 );
  splineY.resize( 0 );
  splineZ.resize( 0 );

  dsplineX.resize( 0 );
  dsplineY.resize( 0 );
  dsplineZ.resize( 0 );
}

Vector3D& 
SplineWarpXform::GetDeformedControlPointPosition
( Vector3D& v, const int x, const int y, const int z) 
  const 
{
  // Create a pointer to the front-lower-left corner of the c.p.g. cell.
  const Types::Coordinate* coeff = m_Parameters + 3 * ( (x-1) + Dims[0] * ((y-1) + Dims[1] * (z-1)) );  
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
    v.XYZ[dim] = mm;
    ++coeff;
    }
  
  return v;
}

void
SplineWarpXform
::GetTransformedGridNonVirtual 
( Vector3D& v, const int idxX, const int idxY, const int idxZ ) const
{
  const Types::Coordinate* coeff = this->m_Parameters + gX[idxX] + gY[idxY] + gZ[idxZ];
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
	  coeff_ll += nextJ;
	  }	
	mm += spZ[m] * ll;
	coeff_mm += nextK;
	}
      v.XYZ[ dim ] = mm;
      ++coeff;
    }
}

void 
SplineWarpXform::GetTransformedGridSequenceNonVirtual
( Vector3D *const vIn, const int numPoints, const int idxX, const int idxY, const int idxZ ) 
  const
{
  Vector3D *v = vIn;
  const Types::Coordinate* coeff = m_Parameters + gX[idxX] + gY[idxY] + gZ[idxZ];
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
  const int numberOfCells = (gX[idxX + numPoints - 1] - gX[idxX]) / nextI + 4;
  
  // pre-compute the contributions of all control points in y- and z-direction
  // along the way
  Types::Coordinate phiComp;
  std::vector<Types::Coordinate> phiHat( 3*numberOfCells );

  const int *gpo;
  int phiIdx = 0;
  for ( int cell = 0; cell < numberOfCells; ++cell, coeff += nextI ) 
    {
    gpo = &GridPointOffset[0];
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

void
SplineWarpXform::ApplyToAll
( CoordinateVector& v, BitVector& valid, const bool inverse, 
 const Types::Coordinate epsilon, const int* gridDims )
{
  const size_t numberOfPoints = v.Dim / 3;
  Types::Coordinate* p = v.Elements;
  
  if ( inverse ) 
    {
    if ( gridDims ) 
      {
      size_t offset = 0;
      Vector3D vv, u;
      for ( int cz = 0; cz < gridDims[2]; ++cz ) 
	{
	for ( int cy = 0; cy < gridDims[1]; ++cy ) 
	  {
	  bool success = false;
	  for ( int cx = 0; cx < gridDims[0]; ++cx, ++offset ) 
	    {
	    if ( !valid[offset] ) continue;
	    vv.Set( p );
	    if ( cx && success )
	      // use previous (successful) inverse as starting point inside row
	      success = this->ApplyInverseInPlaceWithInitial( vv, u, epsilon );
	    else
	      success = this->ApplyInverseInPlace( vv, epsilon );
	    
	    if ( success ) 
	      {
	      u = vv;
	      for ( int dim = 0; dim < 3; ++dim ) p[dim] = vv.XYZ[dim];
	      } 
	    else
	      {
	      valid.Set( offset, false );
	      }
	    }
	  }
	}
      } 
    else
      {
      for ( size_t i = 0; i < numberOfPoints; ++i, p += 3 ) 
	{
	if ( valid[i] ) 
	  {
	  Vector3D vv( p );
	  if ( this->ApplyInverseInPlace( vv, epsilon ) )
	    for ( int dim = 0; dim < 3; ++dim ) p[dim] = vv.XYZ[dim];
	  else
	    valid.Reset( i );
	  }
	}
      }
    } 
  else
    {
    for ( size_t i = 0; i < numberOfPoints; ++i, p += 3 ) 
      {
      if ( valid[i] ) 
	{
	Vector3D vv( p );
	this->ApplyInPlaceNonVirtual( vv );
	for ( int dim = 0; dim < 3; ++dim ) p[dim] = vv.XYZ[dim];
	}
      }
    }
}

Types::Coordinate
SplineWarpXform
::GetGridEnergy() const
{
  double Energy = 0;

  const Types::Coordinate* coeff = m_Parameters + nextI + nextJ + nextK;
  for ( int z = 1; z<Dims[2]-1; ++z, coeff+=2*nextJ )
    for ( int y = 1; y<Dims[1]-1; ++y, coeff+=2*nextI )
      for ( int x = 1; x<Dims[0]-1; ++x, coeff+=nextI ) 
	{
	Energy += this->GetGridEnergy( coeff );
	}
  
  return Energy / (( Dims[0] - 2 ) * ( Dims[1] - 2 ) * ( Dims[2] - 2 ));
}

Types::Coordinate
SplineWarpXform
::GetGridEnergy( const Vector3D& v ) const
{
  Types::Coordinate r[3], f[3];
  int grid[3];
  
  for ( int dim = 0; dim<3; ++dim ) 
    {
    r[dim] = InverseSpacing[dim] * v.XYZ[dim];
    grid[dim] = std::min( static_cast<int>( r[dim] ), Dims[dim]-4 );
    f[dim] = r[dim] - grid[dim];
    assert( (f[dim] >= 0.0) && (f[dim] <= 1.0) );
    }
  
  const Types::Coordinate* coeff = this->m_Parameters + 3 * ( grid[0] + Dims[0] * (grid[1] + Dims[1] * grid[2]) );

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
    MathUtil::Square( InverseSpacing[0] ) *
    ( J[0][0] * J[0][0] + J[0][1] * J[0][1] + J[0][2] * J[0][2] ) +
    MathUtil::Square( InverseSpacing[1] ) *
    ( J[1][0] * J[1][0] + J[1][1] * J[1][1] + J[1][2] * J[1][2] ) +
    MathUtil::Square( InverseSpacing[2] ) *
    ( J[2][0] * J[2][0] + J[2][1] * J[2][1] + J[2][2] * J[2][2] ) +
    // two-variable mixed derivatives
    2 * ( InverseSpacing[0] * InverseSpacing[1] *
	  ( K[0][0] * K[0][0] + K[0][1] * K[0][1] + K[0][2] * K[0][2] ) +
	  InverseSpacing[1] * InverseSpacing[2] *
	  ( K[1][0] * K[1][0] + K[1][1] * K[1][1] + K[1][2] * K[1][2] ) +
	  InverseSpacing[2] * InverseSpacing[0] *
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
    MathUtil::Square( InverseSpacing[0] ) *
    ( J[0][0] * J[0][0] + J[0][1] * J[0][1] + J[0][2] * J[0][2] ) +
    MathUtil::Square( InverseSpacing[1] ) *
    ( J[1][0] * J[1][0] + J[1][1] * J[1][1] + J[1][2] * J[1][2] ) +
    MathUtil::Square( InverseSpacing[2] ) *
    ( J[2][0] * J[2][0] + J[2][1] * J[2][1] + J[2][2] * J[2][2] ) +
    // two-variable mixed derivatives
    2 * ( InverseSpacing[0] * InverseSpacing[1] *
	  ( K[0][0] * K[0][0] + K[0][1] * K[0][1] + K[0][2] * K[0][2] ) +
	  InverseSpacing[1] * InverseSpacing[2] *
	  ( K[1][0] * K[1][0] + K[1][1] * K[1][1] + K[1][2] * K[1][2] ) +
	  InverseSpacing[2] * InverseSpacing[0] *
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
  const unsigned short x =  ( controlPointIdx %  Dims[0] );
  const unsigned short y = ( (controlPointIdx /  Dims[0]) % Dims[1] );
  const unsigned short z = ( (controlPointIdx /  Dims[0]) / Dims[1] );
  
  const int thisDim = param % nextI;
  const Types::Coordinate* coeff = m_Parameters + param - thisDim;
  
  double ground = 0;

  const int iFrom = std::max<int>( -1, 1-x );
  const int jFrom = std::max<int>( -1, 1-y );
  const int kFrom = std::max<int>( -1, 1-z );

  const int iTo = std::min<int>( 1, Dims[0]-2-x );
  const int jTo = std::min<int>( 1, Dims[1]-2-y );
  const int kTo = std::min<int>( 1, Dims[2]-2-z );

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

  upper /= NumberOfControlPoints;
  lower /= NumberOfControlPoints;
}

Types::Coordinate
SplineWarpXform::GetInverseConsistencyError
( const WarpXform* inverse, const UniformVolume* volume,
  const Rect3D* voi )
  const 
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

  const int dX = 1 + static_cast<int>( this->Spacing[0] / 2 * volume->Delta[AXIS_X] );
  const int dY = 1 + static_cast<int>( this->Spacing[1] / 2 * volume->Delta[AXIS_Y] );
  const int dZ = 1 + static_cast<int>( this->Spacing[2] / 2 * volume->Delta[AXIS_Z] );

  const int startX = pVoi->startX - (pVoi->startX % dX);
  const int startY = pVoi->startY - (pVoi->startY % dY);
  const int startZ = pVoi->startZ - (pVoi->startZ % dZ);

  const size_t length = pVoi->endX - startX;
  Array<Vector3D> vecArray( length );

  for ( int z = startZ; z < pVoi->endZ; z += dZ ) 
    {
    for ( int y = startY; y < pVoi->endY; y += dY ) 
      {
      Vector3D* pVec = &vecArray[0];
      this->GetTransformedGridSequenceNonVirtual( pVec, length, startX, y, z );

      for ( int x = startX; x < pVoi->endX; x += dX, pVec += dX ) 
	{
	if ( inverse->InDomain( *pVec ) ) 
	  {
	  inverse->ApplyInPlace( *pVec );
	  volume->GetGridLocation( v, x, y, z );
	  v -= *pVec;
	  result += v.EuclidNorm();
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
  Types::Coordinate* points = Memory::AllocateArray<Types::Coordinate>(  numberOfParameters  );
  memcpy( points, this->m_Parameters, sizeof( *points ) * numberOfParameters );
  
  AffineXform::SmartPtr xform( this->GetInitialAffineXform()->MakeInverse() );

  if ( includeScale ) 
    {
    xform->SetScales( 1.0, 1.0, 1.0 );
    }

  Types::Coordinate* ptr = points;
  Vector3D u;
  for ( size_t pointIdx = 0; pointIdx < numberOfParameters / 3; ++pointIdx, ptr += 3 ) 
    {
    Vector3D v( ptr );
    
    // undo affine transformation component
    xform->ApplyInPlaceNonVirtual( v );
    
    // copy the result into ouput array
    for ( unsigned int dim = 0; dim < 3; ++dim ) 
      ptr[dim] = v.XYZ[dim];
    }
  
  return points;
}

}
