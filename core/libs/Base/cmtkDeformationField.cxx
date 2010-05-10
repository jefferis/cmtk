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

#include <cmtkDeformationField.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
DeformationField::InitControlPoints( const AffineXform* affineXform )
{
  this->m_ParameterVector->Clear();
  
  if ( affineXform ) 
    {
    Types::Coordinate *ofs = this->m_Parameters;

    Self::SpaceVectorType p;
    p[2] = this->m_Offset[2];
    for ( int z = 0; z < this->m_Dims[2]; ++z, p[2] += this->Spacing[2] ) 
      {
      p[1] = this->m_Offset[1];
      for ( int y = 0; y < this->m_Dims[1]; ++y, p[1] += this->Spacing[1] ) 
	{
	p[0] = this->m_Offset[0];
	for ( int x = 0; x < this->m_Dims[0]; ++x, p[0] += this->Spacing[0], ofs+=3 ) 
	  {
	  Self::SpaceVectorType q( p );
	  affineXform->ApplyInPlace( q );
	  q -= p;

	  ofs[0] = q[0];
	  ofs[1] = q[1];
	  ofs[2] = q[2];
	  }
	}
      }	
    
    affineXform->GetScales( this->InverseAffineScaling );
    this->GlobalScaling = affineXform->GetGlobalScaling();
    } 
  else
    {
    this->InverseAffineScaling[0] = this->InverseAffineScaling[1] = this->InverseAffineScaling[2] = this->GlobalScaling = 1.0;
    }
}

void
DeformationField
::GetTransformedGrid 
( Self::SpaceVectorType& v, const int idxX, const int idxY, const int idxZ ) const
{
  const Types::Coordinate* coeff = this->m_Parameters + nextI * idxX + nextJ * idxY + nextK * idxZ;

  v[0] = this->m_Offset[0] + this->Spacing[0] * idxX + coeff[0];
  v[1] = this->m_Offset[1] + this->Spacing[1] * idxY + coeff[1];
  v[2] = this->m_Offset[2] + this->Spacing[2] * idxZ + coeff[2];
}

void 
DeformationField::GetTransformedGridSequence
( Self::SpaceVectorType *const vIn, const int numPoints, const int idxX, const int idxY, const int idxZ ) 
  const
{
  Self::SpaceVectorType *v = vIn;
  const Types::Coordinate* coeff = this->m_Parameters + 3 * (idxX + nextJ * (idxY + nextK * idxZ ));

  const Types::Coordinate Y = this->m_Offset[1] + this->Spacing[1] * idxY;
  const Types::Coordinate Z = this->m_Offset[2] + this->Spacing[2] * idxZ;

  for ( int n = 0; n < numPoints; ++n, ++v, coeff += 3 )
    {
    v[n][0] = this->m_Offset[0] + this->Spacing[0] * idxX + coeff[0];
    v[n][1] = Y + coeff[1];
    v[n][2] = Z + coeff[2];
    }
}

void
DeformationField::ApplyInPlace
( Self::SpaceVectorType& v ) const
{
  Types::Coordinate r[3], f[3];
  int grid[3];
  
  // Do some precomputations.
  for ( int dim = 0; dim<3; ++dim ) 
    {
    // This is the (real-valued) index of the control point grid cell the
    // given location is in.
    r[dim] = this->InverseSpacing[dim] * ( v[dim] - this->m_Offset[dim] );
    // This is the actual cell index.
    grid[dim] = std::min( static_cast<int>( r[dim] ), this->m_Dims[dim]-2 );
    // And here's the relative position within the cell.
    f[dim] = r[dim] - grid[dim];
    }
  
  // Create a pointer to the front-lower-left corner of the c.p.g. cell.
  Types::Coordinate* coeff = this->m_Parameters + 3 * ( grid[0] + this->m_Dims[0] * (grid[1] + this->m_Dims[1] * grid[2]) );
  
  for ( int dim = 0; dim<3; ++dim ) 
    {
    Types::Coordinate mm = 0;
    Types::Coordinate *coeff_mm = coeff;
    
    // Loop over 4 c.p.g. planes in z-direction.
    for ( int m = 0; m < 2; ++m ) 
      {
      Types::Coordinate ll = 0;
      Types::Coordinate *coeff_ll = coeff_mm;
      
      // Loop over 4 c.p.g. planes in y-direction.
      for ( int l = 0; l < 2; ++l ) 
	{
	Types::Coordinate kk = 0;
	Types::Coordinate *coeff_kk = coeff_ll;
	
	// Loop over 4 c.p.g. planes in x-direction.
	for ( int k = 0; k < 2; ++k, coeff_kk+=3 ) 
	  {
	  kk += ( k ? f[0] : 1-f[0] ) * (*coeff_kk);
	  }
	ll += ( l ? f[1] : 1-f[1] ) * kk;
	coeff_ll += nextJ;
	}	
      mm +=  ( m ? f[2] : 1-f[2] ) * ll;
      coeff_mm += nextK;
      }
    v[dim] += mm;
    ++coeff;
    }
}

} // namespace cmtk
