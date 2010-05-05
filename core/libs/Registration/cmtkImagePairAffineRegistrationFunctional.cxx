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

#include <cmtkImagePairAffineRegistrationFunctional.h>

#include <cmtkImagePairAffineRegistrationFunctionalTemplate.h>
#include <cmtkImagePairSimilarityMeasureCR.h>
#include <cmtkImagePairSimilarityMeasureMSD.h>
#include <cmtkImagePairSimilarityMeasureNCC.h>
#include <cmtkImagePairSimilarityMeasureNMI.h>
#include <cmtkImagePairSimilarityMeasureMI.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairAffineRegistrationFunctional* 
ImagePairAffineRegistrationFunctional
::Create
( const int metric, UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume, 
  const Interpolators::InterpolationEnum interpolation, AffineXform::SmartPtr& affineXform )
{
  switch ( metric ) 
    {
    case 0:
      return new ImagePairAffineRegistrationFunctionalTemplate<ImagePairSimilarityMeasureNMI>( refVolume, fltVolume, interpolation, affineXform );
    case 1:
      return new ImagePairAffineRegistrationFunctionalTemplate<ImagePairSimilarityMeasureMI>( refVolume, fltVolume, interpolation, affineXform );
    case 2:
      return new ImagePairAffineRegistrationFunctionalTemplate<ImagePairSimilarityMeasureCR>( refVolume, fltVolume, interpolation, affineXform );
    case 3:
      return NULL; // masked nmi retired
    case 4:
      return new ImagePairAffineRegistrationFunctionalTemplate<ImagePairSimilarityMeasureMSD>( refVolume, fltVolume, interpolation, affineXform );
    case 5:
      return new ImagePairAffineRegistrationFunctionalTemplate<ImagePairSimilarityMeasureNCC>( refVolume, fltVolume, interpolation, affineXform );
    default:
      break;
    }
  return NULL;
}

} // namespace cmtk
