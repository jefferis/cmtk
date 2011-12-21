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

#include <vector>

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

  /// Sinc image interpolation kernel radius.
  static const int InterpolationKernelRadius = 3; 

  /// Constructor.
  EchoPlanarUnwarpFunctional( UniformVolume::SmartConstPtr& imageFwd /*!< "Forward direction" EPI */, 
			      UniformVolume::SmartConstPtr& imageRev /*!< "Reverse" direction EPI */, 
			      const byte phaseEncodeDirection /*!< Phase encoding direction (image coordinate axis) */ )
    : m_ImageFwd( imageFwd ), m_ImageRev( imageRev ), m_PhaseEncodeDirection( phaseEncodeDirection )
  {
    this->m_Deformation.resize( this->m_ImageFwd->GetNumberOfPixels(), 0.0 );
    this->m_UnwarpImageFwd.resize( this->m_ImageFwd->GetNumberOfPixels() );
    this->m_UnwarpImageRev.resize( this->m_ImageFwd->GetNumberOfPixels() );
  }

  /// Return either first or second corrected image.
  UniformVolume::SmartPtr GetCorrectedImage( const byte idx = 0 )
  {
    return UniformVolume::SmartPtr( this->m_ImageFwd->CloneGrid() );
  }
  
private:
  /// "Forward" phase encoding image.
  UniformVolume::SmartConstPtr m_ImageFwd;

  /// "Reverse" phase encoding image.
  UniformVolume::SmartConstPtr m_ImageRev;

  /// Phase encoding direction.
  byte m_PhaseEncodeDirection;

  /// 1D deformation map along phase encoding direction.
  std::vector<Types::Coordinate> m_Deformation;

  /// Deformed "forward" image.
  std::vector<Types::DataItem> m_UnwarpImageFwd;

  /// Deformed "reverse" image.
  std::vector<Types::DataItem> m_UnwarpImageRev;
  
  /// Compute deformed image.
  void ComputeDeformedImage( const UniformVolume& sourceImage /*!< Undeformed input image.*/, UniformVolume& targetImage /*!< Reference to deformed output image.*/, 
			     int direction /*!< Deformation direction - 1 computes unwarped "forward" image, -1 computed unwarped "reverse" image.*/ );

  /// 1D sinc interpolation
  Types::DataItem Interpolate1D( const UniformVolume& sourceImage /*!< Image to interpolate from */, 
				 const FixedVector<3,int>& baseIdx /*!< Grid base index - this is the grid cell where the interpolation kernel is anchored. */, 
				 const Types::Coordinate relative /*!< Relative position of interpolation location in grid cell, relative to phase encoding direction. */ ) const;

  /** Get partial image deformation Jacobian.
   *\return The forward image Jacobian can be computed from this as Jfwd = 1+Jpartial, the reverse Jacobian as Jrev = 1-Jpartial.
   */
  Types::Coordinate GetPartialJacobian( const FixedVector<3,int>& baseIdx /*!< Grid base index - this is the grid cell where the differential operator stencil is anchored. */ ) const;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkEchoPlanarUnwarpFunctional_h_included_
