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

#ifndef __cmtkTransformChangeFromSpaceAffine_h_included_
#define __cmtkTransformChangeFromSpaceAffine_h_included_

#include <cmtkconfig.h>

#include <cmtkAffineXform.h>
#include <cmtkUniformVolume.h>

#include <string>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Compute affine coordinate transformation in standard space from transformation in natrive reference and floating image coordinate spaces.
class TransformChangeFromSpaceAffine
{
public:
  /// Constructor: compute transformation between new spaces.
  TransformChangeFromSpaceAffine( const AffineXform* xform, //!< Transformation from reference to floating in their current spaces.
			      const UniformVolume* reference, //!< Reference (fixed) image.
			      const std::string& referenceSpaceNew, //!< New space for the reference image.
			      const UniformVolume* floating, //! Floating (moving) image.
			      const std::string& floatingSpaceNew //!< New space for the reference image.
    )
  {
    this->ComputeNewTransformation( xform, reference, referenceSpaceNew, floating, floatingSpaceNew );
  }
  
  /// Simplified constructor: compute transformation between images in new, common space.
  TransformChangeFromSpaceAffine( const AffineXform* xform, //!< Transformation from reference to floating in their current spaces.
			      const UniformVolume* reference, //!< Reference (fixed) image.
			      const UniformVolume* floating, //! Floating (moving) image.
			      const std::string& spaceNew //!< New space for both reference and floating image.
    )
  {
    this->ComputeNewTransformation( xform, reference, spaceNew, floating, spaceNew );
  }
  
  /// Simplified constructor: compute transformation between native spaces.
  TransformChangeFromSpaceAffine( const AffineXform* xform, //!< Transformation from reference to floating in their current spaces.
			      const UniformVolume* reference, //!< Reference (fixed) image.
			      const UniformVolume* floating //! Floating (moving) image.
    )
  {
    this->ComputeNewTransformation( xform, reference, reference->m_MetaInformation[CMTK_META_SPACE_ORIGINAL], floating, floating->m_MetaInformation[CMTK_META_SPACE_ORIGINAL] );
  }
  
  /// Return transformation in native spaces.
  const AffineXform& GetTransformation() const
  {
    return this->m_NewXform;
  }

private:
  /// Transformation between native spaces.
  AffineXform m_NewXform;

  /// Actual function to compute transformation between new spaces.
  void ComputeNewTransformation( const AffineXform* xform, //!< Transformation from reference to floating in their current spaces.
				 const UniformVolume* reference, //!< Reference (fixed) image.
				 const std::string& referenceSpaceNew, //!< New space for the reference image.
				 const UniformVolume* floating, //! Floating (moving) image.
				 const std::string& floatingSpaceNew //!< New space for the reference image.
    );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTransformChangeFromSpaceAffine_h_included_
