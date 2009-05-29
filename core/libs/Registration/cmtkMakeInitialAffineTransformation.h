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

#ifndef __cmtkMakeInitialAffineTransformation_h_included_
#define __cmtkMakeInitialAffineTransformation_h_included_

#include <cmtkconfig.h>

#include <cmtkUniformVolume.h>
#include <cmtkAffineXform.h>

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

  /** Align images based on their direction vectors.
   * The direction vectors are encoded in each volume's "AffineXform" field.
   */
  static AffineXform* AlignDirectionVectors
  ( const UniformVolume& referenceImage, //!< The reference (fixed) image
    const UniformVolume& floatingImage, //!< The floating (moving) image
    const bool centerXform = false //!< If this flag is set, the rotation center of the transformation is set to the center of the reference image.
    );

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
