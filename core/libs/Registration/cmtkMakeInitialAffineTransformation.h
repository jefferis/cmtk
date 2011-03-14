/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkMakeInitialAffineTransformation_h_included_
#define __cmtkMakeInitialAffineTransformation_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkAffineXform.h>

#include <string>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Class for generating initial affine coordinate transformations between two images
 */
class MakeInitialAffineTransformation
{
public:
  /// This class.
  typedef MakeInitialAffineTransformation Self;

  /// Enum that defines all initialization modes supported by this class.
  typedef enum
  {
    /// No initialization. Usually means use identity transformation.
    NONE = 0,
    /// Align centers of fields of view.
    FOV = 1,
    /// Align centers of mass.
    COM = 2,
    /// Align using principal axes.
    PAX = 3,
    /// Align using physical coordinates, ie., image origins and direction vectors.
    PHYS = 4
  } Mode;

  /// Return a name for each initialization mode.
  static const std::string GetModeName( const Self::Mode mode );

  /// Create an initial affine transformation for two images based on a selected mode.
  static AffineXform* Create( const UniformVolume& referenceImage /*!< The reference (fixed) image*/,
			      const UniformVolume& floatingImage /*!< The floating (moving) image*/,
			      const Self::Mode mode /*!< Selected initialization method.*/ );

  /** Align images based on their direction vectors.
   * The direction vectors are encoded in each volume's "AffineXform" field.
   */
  static AffineXform* AlignDirectionVectors( const UniformVolume& referenceImage /*!< The reference (fixed) image*/,
					     const UniformVolume& floatingImage /*!< The floating (moving) image*/,
					     const bool centerXform = false /*!< If this flag is set, the rotation center of the transformation is set to the center of the reference image.*/ );
  
  /** Align images based on fields of view.
   *\return This function returns a transformation with three degrees of freedom for a translation only, which
   * aligns the centers of field of view for the two input images. If a crop region is defined in an image, the
   * crop region center is used, otherwise the bounding box center.
   */
  static AffineXform* AlignFieldsOfView( const UniformVolume& referenceImage, const UniformVolume& floatingImage );

  /** Align images based on center of mass.
   *\return This function returns a transformation with three degrees of freedom for a translation only.
   */
  static AffineXform* AlignCentersOfMass( const UniformVolume& referenceImage, const UniformVolume& floatingImage );

  /** Rigidly align images based on principal axes.
   * This function implies alignment by translation according to the centers of mass.
   */
  static AffineXform* AlignPrincipalAxes( const UniformVolume& referenceImage, const UniformVolume& floatingImage );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMakeInitialAffineTransformation_h_included_
