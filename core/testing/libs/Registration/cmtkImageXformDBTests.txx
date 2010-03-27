/*
//
//  Copyright 2010 SRI International
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

#include <cmtkImageXformDB.h>

// test database creation
int
testImageXformDBCreate()
{
  cmtk::ImageXformDB db( "imagexform.sqlite" );
  return 0;
}

// test open of existing file
int
testImageXformDBOpen()
{
  cmtk::ImageXformDB db( CMTK_DATADIR "/empty.sqlite", true /*readOnly*/ );
  return 0;
}

// test adding an image to a database
int
testImageXformDBAddImage()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  db.AddImage( "image1.nii" );
  db.AddImage( "image2.nii", "image1.nii" );
  db.AddImage( "image3.nii", "image2.nii" );
  db.AddImage( "image4.nii" );

  return 0;
}

// test adding images, then a transformation to a database
int
testImageXformDBAddImageThenXform()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  db.AddImage( "image1.nii" );
  db.AddImage( "image2.nii" );
  db.AddImage( "image4.nii" );

  if ( ! db.AddXform( "xform12", true /*invertible*/, "image1.nii", "image2.nii" ) ||
       ! db.AddXform( "xform21", false /*invertible*/, "image2.nii", "image1.nii" ) ||
       ! db.AddXform( "xform42", false /*invertible*/, "image4.nii", "image2.nii" ) )
    return 1;

  return 0;
}

// test adding images and transformations to a database
int
testImageXformDBAddImageWithXform()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  if ( ! db.AddXform( "xform12", true /*invertible*/, "image1.nii", "image2.nii" ) ||
       ! db.AddXform( "xform21", false /*invertible*/, "image2.nii", "image1.nii" ) ||       
       ! db.AddXform( "xform42", false /*invertible*/, "image4.nii", "image2.nii" ) )
    return 1;

  return 0;
}

