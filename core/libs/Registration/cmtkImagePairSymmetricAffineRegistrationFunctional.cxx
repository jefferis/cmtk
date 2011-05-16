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

#include "cmtkImagePairSymmetricAffineRegistrationFunctional.h"

#include <Registration/cmtkImagePairSymmetricAffineRegistrationFunctionalTemplate.h>

#include <Registration/cmtkImagePairSimilarityMeasureCR.h>
#include <Registration/cmtkImagePairSimilarityMeasureMSD.h>
#include <Registration/cmtkImagePairSimilarityMeasureNCC.h>
#include <Registration/cmtkImagePairSimilarityMeasureNMI.h>
#include <Registration/cmtkImagePairSimilarityMeasureMI.h>

#include <Base/cmtkInterpolator.h>
#include <Base/cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairSymmetricAffineRegistrationFunctional* 
ImagePairSymmetricAffineRegistrationFunctional
::Create( const int metric, 
	  UniformVolume::SmartPtr& refVolume,  
	  UniformVolume::SmartPtr& fltVolume,
	  const Interpolators::InterpolationEnum interpolation,
	  AffineXform::SmartPtr& affineXform )
{
  switch ( metric ) 
    {
    case 0:
      return new ImagePairSymmetricAffineRegistrationFunctionalTemplate<ImagePairSimilarityMeasureNMI>( refVolume, fltVolume, interpolation, affineXform );
    case 1:
      return new ImagePairSymmetricAffineRegistrationFunctionalTemplate<ImagePairSimilarityMeasureMI>( refVolume, fltVolume, interpolation, affineXform );
    case 2:
      return new ImagePairSymmetricAffineRegistrationFunctionalTemplate<ImagePairSimilarityMeasureCR>( refVolume, fltVolume, interpolation, affineXform );
    case 3:
      return NULL; // masked NMI retired
    case 4:
      return new ImagePairSymmetricAffineRegistrationFunctionalTemplate<ImagePairSimilarityMeasureMSD>( refVolume, fltVolume, interpolation, affineXform );
    case 5:
      return new ImagePairSymmetricAffineRegistrationFunctionalTemplate<ImagePairSimilarityMeasureNCC>( refVolume, fltVolume, interpolation, affineXform );
    default:
      return NULL;
    }

  return NULL;
}

} // namespace cmtk
