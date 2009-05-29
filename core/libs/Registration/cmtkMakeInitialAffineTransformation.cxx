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

#include <cmtkMakeInitialAffineTransformation.h>

#include <cmtkAffineXform.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

AffineXform* 
MakeInitialAffineTransformation
::AlignDirectionVectors( const UniformVolume& referenceImage, const UniformVolume& floatingImage, const bool centerXform )
{
  if ( referenceImage.m_MetaInformation[CMTK_META_SPACE] != floatingImage.m_MetaInformation[CMTK_META_SPACE] )
    {
    StdErr << "ERROR: coordinate spaces '" << referenceImage.m_MetaInformation[CMTK_META_SPACE]
	      << "' and '" << floatingImage.m_MetaInformation[CMTK_META_SPACE] << "' do not match.\n";
    return NULL;
    }
  
  if ( referenceImage.m_MetaInformation[CMTK_META_EXTERNAL_SPACE_ID] != floatingImage.m_MetaInformation[CMTK_META_EXTERNAL_SPACE_ID] )
    {
    StdErr << "ERROR: semantic coordinate spaces '" << referenceImage.m_MetaInformation[CMTK_META_EXTERNAL_SPACE_ID]
	      << "' and '" << floatingImage.m_MetaInformation[CMTK_META_EXTERNAL_SPACE_ID] << "' do not match.\n";
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
    xform->ChangeCenter( center.XYZ );
    }

  return xform;
}

AffineXform* 
MakeInitialAffineTransformation
::AlignCentersOfMass( const UniformVolume& referenceImage, const UniformVolume& floatingImage )
{
  AffineXform* xform = new AffineXform;
  
  const Vector3D translation = floatingImage.GetCenterOfMass() - referenceImage.GetCenterOfMass();
  xform->SetXlate( translation.XYZ );

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

  // detect pairs of unnecessary 180 deg rotations and fix
  Matrix3x3<Types::Coordinate> product = pAxesRef * pAxesFlt;
  
  int flipCount = 0;
  for ( int i = 0; i < 3; ++i )
    if ( product[i][i] < 0 )
      ++flipCount;
  
  if ( flipCount > 1 )
    {
    int flipA, flipB;
    if ( product[0][0] < product[1][1] )
      {
      flipA = 0;
      if ( product[2][2] < product[1][1] )
	flipB = 2;
      else
	flipB = 1;
      }
    else
      {
      flipA = 1;
      if ( product[0][0] < product[2][2] )
	flipB = 0;
      else
	flipB = 2;
      }
    
    for ( int i = 0; i < 3; ++i )
      {
      pAxesFlt[flipA][i] *= -1;
      pAxesFlt[flipB][i] *= -1;
      }
    }

  // Now compute transformation
  pAxesRef.Invert3x3();
  Matrix3x3<Types::Coordinate> xform3x3 = (pAxesRef * pAxesFlt);

  Vector3D xlation = centerOfMassRef;
  xform3x3.Multiply( xlation.XYZ );
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
  xform->ChangeCenter( centerOfMassRef.XYZ );

  return xform;
}

} // namespace cmtk
