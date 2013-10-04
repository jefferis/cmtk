/*
//
//  Copyright 2004-2013 SRI International
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#include "cmtkVolumeFromFile.h"

#include <Base/cmtkTypedArray.h>

#include <IO/cmtkDICOM.h>

namespace
cmtk
{

const UniformVolume::SmartPtr
VolumeFromFile::ReadDICOM( const std::string& path )
{
  DICOM dicom;
  try
    {
    dicom.InitFromFile( path );
    }
  catch ( const cmtk::Exception& ex )
    {
    StdErr << "ERROR: " << ex.what() << "\n";
    return UniformVolume::SmartPtr( NULL );
    }

  FixedVector<3,int> dims = dicom.GetDims();
  FixedVector<3,double> pixelSize = dicom.GetPixelSize();
  
  const unsigned long totalImageSizePixels = dims[0] * dims[1] * dims[2];
  
  TypedArray::SmartPtr pixelDataArray = dicom.GetPixelDataArray( totalImageSizePixels );
  
  UniformVolume::CoordinateVectorType imageOrigin = dicom.GetImageOrigin();
  FixedArray< 2, FixedVector<3,double> > imageOrientation = dicom.GetImageOrientation();
  
  // without further information, we "guess" the image normal vector
  UniformVolume::CoordinateVectorType sliceNormal = dicom.DemosaicAndGetNormal( imageOrientation, pixelSize, dims, pixelDataArray, imageOrigin );   
    
  // Construct volume and set the DICOM coordinates
  UniformVolume::SmartPtr volume( new UniformVolume( UniformVolume::IndexType( dims ), pixelSize[0], pixelSize[1], pixelSize[2], pixelDataArray ) );
  volume->SetMetaInfo( META_SPACE, "LPS" );
  volume->SetMetaInfo( META_SPACE_ORIGINAL, "LPS" );

  imageOrientation[0] *= pixelSize[0] / imageOrientation[0].RootSumOfSquares();
  imageOrientation[1] *= pixelSize[1] / imageOrientation[1].RootSumOfSquares();
  sliceNormal *= pixelSize[2] / sliceNormal.RootSumOfSquares();

  const Types::Coordinate directions[3][3] = 
    {
      { imageOrientation[0][0], imageOrientation[0][1], imageOrientation[0][2] },
      { imageOrientation[1][0], imageOrientation[1][1], imageOrientation[1][2] },
      { sliceNormal[0], sliceNormal[1], sliceNormal[2] }
    };
  
  const Matrix3x3<Types::Coordinate> m3( directions );
  Matrix4x4<Types::Coordinate> m4( m3 );
  for ( int i = 0; i < 3; ++i )
    m4[3][i] = imageOrigin[i];

  volume->m_IndexToPhysicalMatrix = m4;
//  const std::string orientationString0 = volume->GetOrientationFromDirections();
  volume->ChangeCoordinateSpace( AnatomicalOrientation::ORIENTATION_STANDARD );

  const std::string orientationString = volume->GetOrientationFromDirections();
  volume->SetMetaInfo( META_SPACE_UNITS_STRING, "mm" ); // seems to be implied in DICOM
  volume->SetMetaInfo( META_IMAGE_ORIENTATION, orientationString );
  volume->SetMetaInfo( META_IMAGE_ORIENTATION_ORIGINAL, orientationString );

  return volume;
}

} // namespace cmtk
