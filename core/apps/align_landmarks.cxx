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

#include <cmtkconfig.h>

#include <vtkLandmarkTransform.h>
#include <vtkMatrix4x4.h>

#include <cmtkVolumeIO.h>
#include <cmtkStudy.h>
#include <cmtkStudyList.h>
#include <cmtkClassStreamStudyList.h>

#include <cmtkAffineXform.h>

#include <cmtkLandmarkListVTK.h>

int main ( const int argc, const char *argv[] )
{
  if ( argc < 3 )
    return 1;

  cmtk::Study::SmartPtr refStudy( cmtk::Study::Read( argv[1] ) );
  cmtk::Study::SmartPtr fltStudy( cmtk::Study::Read( argv[2] ) );

  cmtk::LandmarkListVTK::SmartPtr refLL( new cmtk::LandmarkListVTK( *(refStudy->GetLandmarkList()) ) );
  cmtk::LandmarkListVTK::SmartPtr fltLL( new cmtk::LandmarkListVTK( *(fltStudy->GetLandmarkList()) ) );

  vtkPoints *refPoints, *fltPoints;
  refPoints = refLL->GetMatchedVtkPoints( fltPoints, fltLL );

  vtkLandmarkTransform *transform = vtkLandmarkTransform::New();
  transform->SetSourceLandmarks( fltPoints );
  transform->SetTargetLandmarks( refPoints );

  vtkMatrix4x4* matrix = vtkMatrix4x4::New();
  transform->GetMatrix( matrix );
  matrix->Transpose();

  vtkPoints *xfPoints = vtkPoints::New();
  transform->TransformPoints( fltPoints, xfPoints );

  float sqError = 0;
  cmtk::Vector3D u, v;
  const vtkIdType numPoints = refPoints->GetNumberOfPoints();
  for ( vtkIdType id = 0; id < numPoints; ++id ) 
    {
    refPoints->GetPoint( id, u.begin() );
    xfPoints->GetPoint( id, v.begin() );
    u -= v;
    sqError += u.RootSumOfSquares();
    }
  fprintf( stderr, "Matched %d points. FRE = %f [mm].\n", static_cast<int>( numPoints ), static_cast<float>( sqError / numPoints ) );
  
  cmtk::AffineXform::SmartPtr affineXform( new cmtk::AffineXform( matrix->Element ) );
  cmtk::WarpXform::SmartPtr warpXform( NULL );
  
  cmtk::StudyList::SmartPtr studyList( new cmtk::StudyList );
  studyList->AddStudy( refStudy );
  studyList->AddStudy( fltStudy );
  studyList->AddXform( refStudy, fltStudy, affineXform, warpXform );
  
  cmtk::ClassStreamStudyList::Write( argv[3], studyList);
  
  matrix->Delete();
  transform->Delete();
  
  return 0;
}

