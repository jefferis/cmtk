/*
//
//  Copyright 2011 SRI International
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

#ifndef __cmtkEchoPlanarUnwarpFunctional_h_included_
#define __cmtkEchoPlanarUnwarpFunctional_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Class for computing a const function for unwarping Echo Planar Images.
 *\see D. Holland, J. M. Kuperman, and A. M. Dale, "Efficient correction of inhomogeneous static magnetic field-induced distortion in Echo Planar Imaging," NeuroImage, vol. 50, no. 1, pp. 175-183, 2010.
 */
class EchoPlanarUnwarpFunctional
{
public:
  /// This class.
  typedef EchoPlanarUnwarpFunctional Self;

  /// Constructor.
  EchoPlanarUnwarpFunctional( UniformVolume::SmartConstPtr& image0, UniformVolume::SmartConstPtr& image1, const byte phaseEncodeDirection )
    : m_Image0( image0 ), m_Image1( image1 ), m_PhaseEncodeDirection( phaseEncodeDirection )
  {
    this->m_Deformation = image0->CloneGrid();
    this->m_Deformation->CreateDataArray( TYPE_COORDINATE );
    
    this->m_UnwarpImage0 = image0->CloneGrid();
    this->m_UnwarpImage0->CreateDataArray( TYPE_FLOAT );

    this->m_UnwarpImage1 = image1->CloneGrid();
    this->m_UnwarpImage1->CreateDataArray( TYPE_FLOAT );
  }

  /// Return either first or second corrected image.
  const UniformVolume::SmartPtr& GetCorrectedImage( const byte idx = 0 )
  {
    if ( idx == 0 )
      return this->m_UnwarpImage0;
    else
      return this->m_UnwarpImage1;
  }
  
private:
  /// Phase encoding direction.
  byte m_PhaseEncodeDirection;

  /// First image.
  UniformVolume::SmartConstPtr m_Image0;

  /// Second image.
  UniformVolume::SmartConstPtr m_Image1;

  /// Deformation map.
  UniformVolume::SmartPtr m_Deformation;

  /// Deformed first image.
  UniformVolume::SmartPtr m_UnwarpImage0;

  /// Deformed second image.
  UniformVolume::SmartPtr m_UnwarpImage1;

};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkEchoPlanarUnwarpFunctional_h_included_
