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

#include <cmtkUniformVolume.h>

void
cmtk::UniformVolume
::GetGradientAt( Vector3D& g, const int i, const int j, const int k )
{
  g.XYZ[0] = (this->GetDataAt( i+1, j, k ) - this->GetDataAt( i-1, j, k )) / (2*this->m_Delta[0]);
  g.XYZ[1] = (this->GetDataAt( i, j+1, k ) - this->GetDataAt( i, j-1, k )) / (2*this->m_Delta[1]);
  g.XYZ[2] = (this->GetDataAt( i, j, k+1 ) - this->GetDataAt( i, j, k-1 )) / (2*this->m_Delta[2]);
}

void
cmtk::UniformVolume
::GetHessianAt( Matrix3x3<Types::DataItem>& H, const int i, const int j, const int k )
{
// implementation following central differences formulas from http://www.technion.ac.il/docs/sas/ormp/chap5/sect28.htm

  const Types::DataItem central = 30 * this->GetDataAt( i, j, k );
  H[0][0] = (-this->GetDataAt( i+2, j, k ) + 16 * this->GetDataAt( i+1, j, k ) - central + 16 * this->GetDataAt( i-1, j, k ) - this->GetDataAt( i-2, j, k ))
    / ( 12 * this->m_Delta[0] * this->m_Delta[0] );
  H[1][1] = (-this->GetDataAt( i, j+2, k ) + 16 * this->GetDataAt( i, j+1, k ) - central + 16 * this->GetDataAt( i, j-1, k ) - this->GetDataAt( i, j-2, k ))
    / ( 12 * this->m_Delta[1] * this->m_Delta[1] );
  H[2][2] = (-this->GetDataAt( i, j, k+2 ) + 16 * this->GetDataAt( i, j, k+1 ) - central + 16 * this->GetDataAt( i, j, k-1 ) - this->GetDataAt( i, j, k-2 ))
    / ( 12 * this->m_Delta[2] * this->m_Delta[2] );
  
  H[0][1] = H[1][0] = (this->GetDataAt( i+1, j+1, k ) - this->GetDataAt( i+1, j-1, k ) - this->GetDataAt( i-1, j+1, k ) + this->GetDataAt( i-1, j-1, k )) / ( 4 * this->m_Delta[0] * this->m_Delta[1] );
  H[0][2] = H[2][0] = (this->GetDataAt( i+1, j, k+1 ) - this->GetDataAt( i+1, j, k-1 ) - this->GetDataAt( i-1, j, k+1 ) + this->GetDataAt( i-1, j, k-1 )) / ( 4 * this->m_Delta[0] * this->m_Delta[2] );
  H[1][2] = H[2][1] = (this->GetDataAt( i, j+1, k+1 ) - this->GetDataAt( i, j+1, k-1 ) - this->GetDataAt( i, j-1, k+1 ) + this->GetDataAt( i, j-1, k-1 )) / ( 4 * this->m_Delta[1] * this->m_Delta[2] );
}
