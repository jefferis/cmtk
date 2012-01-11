/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include <Registration/cmtkImagePairAffineRegistrationFunctional.h>

#include <Registration/cmtkImagePairAffineRegistrationFunctionalTemplate.h>
#include <Registration/cmtkImagePairSimilarityMeasureCR.h>
#include <Registration/cmtkImagePairSimilarityMeasureMSD.h>
#include <Registration/cmtkImagePairSimilarityMeasureRMS.h>
#include <Registration/cmtkImagePairSimilarityMeasureNCC.h>
#include <Registration/cmtkImagePairSimilarityMeasureNMI.h>
#include <Registration/cmtkImagePairSimilarityMeasureMI.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

Types::Coordinate
ImagePairAffineRegistrationFunctional::GetParamStep( const size_t idx, const Types::Coordinate mmStep ) const 
{
  switch ( this->m_RestrictToInPlane )
    {
    case 0:
      switch ( idx )
	{
	case 0: //xlate x
	case 4: //rot y
	case 5: //rot z
	case 6: //scale x
	case 9: //shear xy
	case 10: //shear xz
	  return 0.0;
	default:
	  break;
	}
      break;
    case 1:
      switch ( idx )
	{
	case 1: //xlate y
	case 3: //rot x
	case 5: //rot z
	case 7: //scale y
	case 9: //shear xy
	case 11: //shear yz
	  return 0.0;
	default:
	  break;
	}
      break;
    case 2:
      switch ( idx )
	{
	case 2: //xlate z
	case 3: //rot x
	case 4: //rot y
	case 8: //scale z
	case 10: //shear xz
	case 11: //shear yz
	  return 0.0;
	default:
	  break;
	}
      break;
    default:
      break;
    }
  
  return this->m_AffineXform->GetParamStep( idx, this->m_FloatingSize, mmStep );  
}

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
      return new ImagePairAffineRegistrationFunctionalTemplate<ImagePairSimilarityMeasureRMS>( refVolume, fltVolume, interpolation, affineXform );
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
