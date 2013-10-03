/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013 SRI International
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

#include "cmtkTransformChangeToSpaceAffine.h"

#include <System/cmtkExitException.h>

cmtk::TransformChangeToSpaceAffine
::TransformChangeToSpaceAffine
( const AffineXform& xform, const UniformVolume& reference, const UniformVolume& floating, const char* forceSpace )
{
// adapt transformation to Slicer's image coordinate systems as defined in the Nrrd files we probably read
  UniformVolume::SmartPtr refVolumeOriginalSpace( reference.CloneGrid() );
  UniformVolume::SmartPtr fltVolumeOriginalSpace( floating.CloneGrid() );
  
  // first bring volumes back into their native coordinate space, or into forced space if one is provided.
  if ( forceSpace )
    {
    refVolumeOriginalSpace->ChangeCoordinateSpace( forceSpace );
    fltVolumeOriginalSpace->ChangeCoordinateSpace( forceSpace );
    }
  else
    {
    refVolumeOriginalSpace->ChangeCoordinateSpace( reference.GetMetaInfo( META_SPACE_ORIGINAL ) );
    fltVolumeOriginalSpace->ChangeCoordinateSpace( floating.GetMetaInfo( META_SPACE_ORIGINAL ) );
    }
  
  // now determine image-to-physical transformations and concatenate these.
  const AffineXform::MatrixType refMatrix = refVolumeOriginalSpace->GetImageToPhysicalMatrix();
  const AffineXform::MatrixType fltMatrix = fltVolumeOriginalSpace->GetImageToPhysicalMatrix();

  try
    {
    const AffineXform::MatrixType concatMatrix = (refMatrix.GetInverse() * xform.Matrix) * fltMatrix;
    
    // create output transformation and write
    this->m_NewXform.SetMatrix( concatMatrix );
    }
  catch ( const AffineXform::MatrixType::SingularMatrixException& ex )
    {
    StdErr << "ERROR: singular matrix in TransformChangeToSpaceAffine constructor.\n";
    throw ExitException( 1 );
    }
}
