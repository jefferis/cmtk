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

#include "cmtkMakeInitialAffineTransformation.h"

#include <Base/cmtkAffineXform.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

const std::string 
MakeInitialAffineTransformation
::GetModeName( const Self::Mode mode )
{
  switch ( mode )
    {
    case Self::NONE:
    default:
      return std::string( "none" );
    case Self::FOV:
      return std::string( "FieldsOfView" );
    case Self::COM:
      return std::string( "CentersOfMass" );
    case Self::PAX:
      return std::string( "PrincipalAxes" );
    case Self::PHYS:
      return std::string( "PhysicalCoordinates" );
    }
  return std::string( "unknown" );
}

AffineXform* 
MakeInitialAffineTransformation
::Create( const UniformVolume& referenceImage, const UniformVolume& floatingImage, const Self::Mode mode )
{
  switch ( mode )
    {
    case Self::NONE:
    default:
      return new AffineXform;
    case Self::FOV:
      return Self::AlignFieldsOfView( referenceImage, floatingImage );
    case Self::COM:
      return Self::AlignCentersOfMass( referenceImage, floatingImage );
    case Self::PAX:
      return Self::AlignPrincipalAxes( referenceImage, floatingImage );
    case Self::PHYS:
      return Self::AlignDirectionVectors( referenceImage, floatingImage );
    }
  return new AffineXform;
}

AffineXform* 
MakeInitialAffineTransformation
::AlignDirectionVectors( const UniformVolume& referenceImage, const UniformVolume& floatingImage, const bool centerXform )
{
  if ( referenceImage.GetMetaInfo( META_SPACE ) != floatingImage.GetMetaInfo( META_SPACE ) )
    {
    StdErr << "ERROR: coordinate spaces '" << referenceImage.GetMetaInfo( META_SPACE )
	   << "' and '" << floatingImage.GetMetaInfo( META_SPACE ) << "' do not match.\n";
    return NULL;
    }
  
  if ( referenceImage.GetMetaInfo( META_EXTERNAL_SPACE_ID ) != floatingImage.GetMetaInfo( META_EXTERNAL_SPACE_ID ) )
    {
    StdErr << "ERROR: semantic coordinate spaces '" << referenceImage.GetMetaInfo( META_EXTERNAL_SPACE_ID )
	   << "' and '" << floatingImage.GetMetaInfo( META_EXTERNAL_SPACE_ID ) << "' do not match.\n";
    return NULL;
    }
  
  const AffineXform::MatrixType refMatrix = referenceImage.GetImageToPhysicalMatrix();
  AffineXform referenceXform( refMatrix );

  const AffineXform::MatrixType fltMatrix = floatingImage.GetImageToPhysicalMatrix();
  AffineXform floatingXform( fltMatrix );
  
  AffineXform* xform = new AffineXform( referenceXform );
  xform->Concat( *floatingXform.GetInverse() );

  if ( centerXform )
    {
    const Vector3D center = referenceImage.GetCenterCropRegion();
    xform->ChangeCenter( center );
    }

  return xform;
}

AffineXform* 
MakeInitialAffineTransformation
::AlignFieldsOfView( const UniformVolume& referenceImage, const UniformVolume& floatingImage )
{
  AffineXform* xform = new AffineXform;
  
  const Vector3D translation = floatingImage.GetCenterCropRegion() - referenceImage.GetCenterCropRegion();
  xform->SetXlate( translation.begin() );
  
  return xform;
}

AffineXform* 
MakeInitialAffineTransformation
::AlignCentersOfMass( const UniformVolume& referenceImage, const UniformVolume& floatingImage )
{
  AffineXform* xform = new AffineXform;
  
  const Vector3D translation = floatingImage.GetCenterOfMass() - referenceImage.GetCenterOfMass();
  xform->SetXlate( translation.begin() );

  return xform;
}

AffineXform* 
MakeInitialAffineTransformation
::AlignPrincipalAxes( const UniformVolume& referenceImage, const UniformVolume& floatingImage )
{
  // get principal axes
  Matrix3x3<Types::Coordinate> pAxesRef, pAxesFlt;
  Vector3D centerOfMassRef, centerOfMassFlt;

  referenceImage.GetPrincipalAxes( pAxesRef, centerOfMassRef );
  floatingImage.GetPrincipalAxes( pAxesFlt, centerOfMassFlt );

  pAxesRef = pAxesRef.GetTranspose();
  pAxesFlt = pAxesFlt.GetTranspose();

  // Now compute transformation
  pAxesRef.Invert3x3();
  Matrix3x3<Types::Coordinate> xform3x3 = (pAxesRef * pAxesFlt);

  Vector3D xlation = centerOfMassRef;
  xform3x3.Multiply( xlation );
  xlation = centerOfMassFlt - xlation;
  
  // Assign xform3x3 as a submatrix of a 4x4
  Matrix4x4<Types::Coordinate> xform4x4 = xform3x3;
  
  // Turn xform4x4 into homogenized matrix with xlation as the 4th column.
  for ( int i = 0; i < 3; i++ )
    {
    xform4x4[3][i] = xlation[i];
    xform4x4[i][3] = 0;
    }
  xform4x4[3][3] = 1;
  
  AffineXform* xform = new AffineXform( xform4x4 );
  xform->ChangeCenter( centerOfMassRef );

  Types::Coordinate* angles = xform->RetAngles();
  for ( int i = 0; i < 3; ++i )
    {
    if ( angles[i] > 90 )
      angles[i] -= 180;
    else if ( angles[i] < -90 )
      angles[i] += 180;
    }
  xform->SetAngles( angles );
  
  return xform;
}

} // namespace cmtk
