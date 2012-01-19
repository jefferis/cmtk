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

#include <Registration/cmtkGroupwiseRegistrationFunctionalAffineInitializer.h>

#include <Base/cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

void
GroupwiseRegistrationFunctionalAffineInitializer::InitializeXforms
( GroupwiseRegistrationFunctionalBase& functional, const bool alignCenters, const bool alignCenterOfMass, const bool initScales, const bool centerInTemplateFOV )
{
  const size_t numberOfImages = functional.m_ImageVector.size();

  const Vector3D centerTemplate = functional.m_TemplateGrid->GetCenterCropRegion();
  
  std::vector<Vector3D> centers( numberOfImages );
  std::vector<Vector3D> firstOrderMoments;
  if ( initScales )
    firstOrderMoments.resize( numberOfImages );
  functional.m_XformVector.resize( numberOfImages );

  Vector3D centerAverage;
  // first get all image centers (either FOV or center of mass)
  for ( size_t imageIdx = 0; imageIdx < numberOfImages; ++imageIdx )
    {
    if ( alignCenters )
      {
      if ( alignCenterOfMass )
	{
	if ( initScales )
	  {
	  centers[imageIdx] = functional.m_ImageVector[imageIdx]->GetCenterOfMass( firstOrderMoments[imageIdx] );
	  }
	else
	  {
	  centers[imageIdx] = functional.m_ImageVector[imageIdx]->GetCenterOfMass();
	  }
	}
      else
	{
	centers[imageIdx] = functional.m_ImageVector[imageIdx]->GetCenter();
	}
      }
    }
  
  if ( centerInTemplateFOV )
    {
    // if so desired, center in template image FOV
    centerAverage = centerTemplate;
    }
  else
    {
    // if not aligning with template FOV, compute centrioid of input image centers.
    std::fill( centerAverage.begin(), centerAverage.end(), 0 );
    for ( size_t imageIdx = 0; imageIdx < numberOfImages; ++imageIdx )
      {
      centerAverage += centers[imageIdx];
      }
    
    // compute average of all image centers
    centerAverage *= (1.0 / numberOfImages);
    }
  
  // now make sure every image gests shifted so their center align, and the overall shift is zero
  for ( size_t imageIdx = 0; imageIdx < numberOfImages; ++imageIdx )
    {
    AffineXform::SmartPtr xform( new AffineXform );
    xform->SetUseLogScaleFactors( true );
    xform->SetCenter( centerTemplate.begin() );
 
    const Vector3D delta( centers[imageIdx] - centerAverage );
    
    xform->SetXlate( delta.begin() );
    functional.m_XformVector[imageIdx] = xform;
    }
  
  // convert first order moments to scale with average log factor 0
  if ( initScales )
    {
    UniformVolume::CoordinateVectorType avgScales( UniformVolume::CoordinateVectorType::Init( 0 ) );
    UniformVolume::CoordinateVectorType fom0( firstOrderMoments[0] );
    for ( size_t imageIdx = 0; imageIdx < numberOfImages; ++imageIdx )
      {
      for ( int dim = 0; dim < 3; ++dim )
	firstOrderMoments[imageIdx][dim] = log( firstOrderMoments[imageIdx][dim] / fom0[dim] );
      avgScales += firstOrderMoments[imageIdx];
      }
    avgScales *= ( 1.0 /  numberOfImages );
    for ( size_t imageIdx = 0; imageIdx < numberOfImages; ++imageIdx )
      {
      firstOrderMoments[imageIdx] -= avgScales;
      AffineXform::SmartPtr::DynamicCastFrom( functional.m_XformVector[imageIdx] )->SetScales( firstOrderMoments[imageIdx].begin() );
      }
    }
}

} // namespace cmtk
