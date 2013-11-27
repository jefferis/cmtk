/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include "cmtkCompatibilityMatrix4x4.h"

template<class T>
cmtk::CompatibilityMatrix4x4<T>::CompatibilityMatrix4x4( const CoordinateVector& dofs, const bool logScaleFactors )
{
  const Units::Radians alpha = Units::Degrees( dofs[3] );
  const Units::Radians theta = Units::Degrees( dofs[4] );
  const Units::Radians   phi = Units::Degrees( dofs[5] );

  const double cos0 = MathUtil::Cos(alpha), sin0 = MathUtil::Sin(alpha);
  const double cos1 = MathUtil::Cos(theta), sin1 = MathUtil::Sin(theta);
  const double cos2 = MathUtil::Cos(  phi), sin2 = MathUtil::Sin(  phi);

  const double sin0xsin1 = sin0 * sin1;
  const double cos0xsin1 = cos0 * sin1;

  const double scaleX = (logScaleFactors) ? exp( dofs[6] ) : dofs[6];
  const double scaleY = (logScaleFactors) ? exp( dofs[7] ) : dofs[7];
  const double scaleZ = (logScaleFactors) ? exp( dofs[8] ) : dofs[8];

  this->m_Matrix[0][0] = static_cast<T>( cos1*cos2 * scaleX );
  this->m_Matrix[0][1] = static_cast<T>( -cos1*sin2 * scaleX );                     
  this->m_Matrix[0][2] = static_cast<T>( -sin1 * scaleX );
  this->m_Matrix[0][3] = static_cast<T>( 0 );
  this->m_Matrix[1][0] = static_cast<T>(  (sin0xsin1*cos2 + cos0*sin2) * scaleY );
  this->m_Matrix[1][1] = static_cast<T>( (-sin0xsin1*sin2 + cos0*cos2) * scaleY ); 
  this->m_Matrix[1][2] = static_cast<T>(  sin0*cos1 * scaleY );
  this->m_Matrix[1][3] = static_cast<T>( 0 );
  this->m_Matrix[2][0] = static_cast<T>(  (cos0xsin1*cos2 - sin0*sin2) * scaleZ );
  this->m_Matrix[2][1] = static_cast<T>( (-cos0xsin1*sin2 - sin0*cos2) * scaleZ );
  this->m_Matrix[2][2] = static_cast<T>(  cos0*cos1 * scaleZ );
  this->m_Matrix[2][3] = static_cast<T>( 0 );

  this->m_Matrix[3][0] = this->m_Matrix[3][1] = this->m_Matrix[3][2] = static_cast<T>( 0 );
  this->m_Matrix[3][3] = static_cast<T>( 1.0 );

  // generate shears
  for ( int i = 2; i >= 0; --i )
    { 
    Superclass shear = Superclass::Identity();
    shear[i/2][(i/2)+(i%2)+1] = dofs[9+i];
    *this *= shear;
    }
  
  // transform rotation center
  const Types::Coordinate cM[3] = 
    {
      dofs[12]*this->m_Matrix[0][0] + dofs[13]*this->m_Matrix[1][0] + dofs[14]*this->m_Matrix[2][0],
      dofs[12]*this->m_Matrix[0][1] + dofs[13]*this->m_Matrix[1][1] + dofs[14]*this->m_Matrix[2][1],
      dofs[12]*this->m_Matrix[0][2] + dofs[13]*this->m_Matrix[1][2] + dofs[14]*this->m_Matrix[2][2]
    };
  
  // set translations
  this->m_Matrix[3][0] = dofs[0] - cM[0] + dofs[12];
  this->m_Matrix[3][1] = dofs[1] - cM[1] + dofs[13];
  this->m_Matrix[3][2] = dofs[2] - cM[2] + dofs[14];
}

template class cmtk::CompatibilityMatrix4x4<cmtk::Types::Coordinate>;
