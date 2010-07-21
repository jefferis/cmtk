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

#include "cmtkImagePairSymmetricNonrigidRegistrationFunctional.h"

#include "Registration/cmtkImagePairSymmetricNonrigidRegistrationFunctionalTemplate.h"
#include "Registration/cmtkImagePairSimilarityMeasureCR.h"
#include "Registration/cmtkImagePairSimilarityMeasureMSD.h"
#include "Registration/cmtkImagePairSimilarityMeasureNCC.h"
#include "Registration/cmtkImagePairSimilarityMeasureNMI.h"
#include "Registration/cmtkImagePairSimilarityMeasureMI.h"

#include "Base/cmtkInterpolator.h"
#include "Base/cmtkUniformVolume.h"
#include "Base/cmtkSplineWarpXform.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class VM, class W>
void
ImagePairSymmetricNonrigidRegistrationFunctionalTemplate<VM,W>::SetWarpXform
( SplineWarpXform::SmartPtr& warpFwd, SplineWarpXform::SmartPtr& warpBwd ) 
{
  this->FwdFunctional.SetWarpXform( warpFwd );
  this->FwdFunctional.SetInverseTransformation( warpBwd );
  
  this->BwdFunctional.SetWarpXform( warpBwd );
  this->BwdFunctional.SetInverseTransformation( warpFwd );
}

template<class VM, class W>
typename ImagePairSymmetricNonrigidRegistrationFunctionalTemplate<VM,W>::ReturnType
ImagePairSymmetricNonrigidRegistrationFunctionalTemplate<VM,W>
::EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step )
{
  CoordinateVector vFwd( this->FwdFunctional.ParamVectorDim(), v.Elements, false /*freeElements*/ );
  CoordinateVector gFwd( this->FwdFunctional.ParamVectorDim(), g.Elements, false /*freeElements*/ );

  CoordinateVector vBwd( this->BwdFunctional.ParamVectorDim(), v.Elements+this->FwdFunctional.ParamVectorDim(), false /*freeElements*/ );
  CoordinateVector gBwd( this->BwdFunctional.ParamVectorDim(), g.Elements+this->FwdFunctional.ParamVectorDim(), false /*freeElements*/ );

  const typename Self::ReturnType result =
    this->FwdFunctional.EvaluateWithGradient( vFwd, gFwd, step ) + this->BwdFunctional.EvaluateWithGradient( vBwd, gBwd, step );
  return result;
}

ImagePairSymmetricNonrigidRegistrationFunctional* 
ImagePairSymmetricNonrigidRegistrationFunctional
::Create( const int metric, 
	  UniformVolume::SmartPtr& refVolume,  
	  UniformVolume::SmartPtr& fltVolume,
	  const Interpolators::InterpolationEnum interpolation )
{
  switch ( metric ) 
    {
    case 0:
      return new ImagePairSymmetricNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureNMI,SplineWarpXform>( refVolume, fltVolume, interpolation );
    case 1:
      return new ImagePairSymmetricNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureMI,SplineWarpXform>( refVolume, fltVolume, interpolation );
    case 2:
      return new ImagePairSymmetricNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureCR,SplineWarpXform>( refVolume, fltVolume, interpolation );
    case 3:
      return NULL; // masked NMI retired
    case 4:
      return new ImagePairSymmetricNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureMSD,SplineWarpXform>( refVolume, fltVolume, interpolation );
    case 5:
      return new ImagePairSymmetricNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureNCC,SplineWarpXform>( refVolume, fltVolume, interpolation );
    default:
      return NULL;
    }

  return NULL;
}

} // namespace cmtk
