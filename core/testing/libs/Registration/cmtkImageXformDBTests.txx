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

  if ( ! db.AddImagePairXform( "xform12", true /*invertible*/, "image1.nii", "image2.nii" ) ||
       ! db.AddImagePairXform( "xform21", false /*invertible*/, "image2.nii", "image1.nii" ) ||
       ! db.AddImagePairXform( "xform42", false /*invertible*/, "image4.nii", "image2.nii" ) )
    return 1;

  return 0;
}

// test adding images and transformations to a database
int
testImageXformDBAddImageWithXform()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  if ( ! db.AddImagePairXform( "xform12", true /*invertible*/, "image1.nii", "image2.nii" ) ||
       ! db.AddImagePairXform( "xform21", false /*invertible*/, "image2.nii", "image1.nii" ) ||       
       ! db.AddImagePairXform( "xform42", false /*invertible*/, "image4.nii", "image2.nii" ) )
    return 1;

  return 0;
}

// test: try adding a transformation between images in the same space
int
testImageXformDBAddXformSameSpace()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  db.AddImage( "image1.nii" );
  db.AddImage( "image2.nii", "image1.nii" );
  if ( db.AddImagePairXform( "xform12", true /*invertible*/, "image1.nii", "image2.nii" ) )
    return 1;
  
  return 0;
}

// test: add transformation, then add its refinement.
int
testImageXformDBAddXformRefined()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  db.AddImage( "image1.nii" );
  db.AddImage( "image2.nii" );
  if ( ! db.AddImagePairXform( "xform12", true /*invertible*/, "image1.nii", "image2.nii" ) ||
       ! db.AddRefinedXform( "xform12refined", true /*invertible*/, "xform12" ) )
    {
    std::cerr << "DB add transformation." << std::cerr;
    return 1;
    }

  if ( (db.FindXformLevel( "xform12" ) != 0) ||
       (db.FindXformLevel( "xform12refined" ) != 1) )
    {
    std::cerr << "DB transformation levels are incorrect." << std::cerr;
    return 1;
    }
  
  return 0;
}

// test getting simple transformations
int
testImageXformDBFindXform()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  if ( ! db.AddImagePairXform( "xform12", true /*invertible*/, "image1.nii", "image2.nii" ) )
    return 1;
  
  std::string xform;
  bool inverse;
  
  if ( ! db.FindXform( "image1.nii", "image2.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::cerr;
    return 1;
    }

  if ( xform != "xform12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::cerr;
    return 1;
    }

  if ( inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::cerr;
    return 1;
    }

  return 0;
}

// test getting transformation between two images in the same space.
int
testImageXformDBFindXformSameSpace()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();

  db.AddImage( "image1.nii" );
  db.AddImage( "image2.nii", "image1.nii" );

  std::string xform;
  bool inverse;
  
  if ( ! db.FindXform( "image1.nii", "image2.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::cerr;
    return 1;
    }

  if ( xform != "" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::cerr;
    return 1;
    }
  
  if ( inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::cerr;
    return 1;
    }

  return 0;
}

// test getting inverse transformations
int
testImageXformDBFindXformInverse()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  if ( ! db.AddImagePairXform( "xform12", true /*invertible*/, "image1.nii", "image2.nii" ) )
    return 1;
  
  std::string xform;
  bool inverse;
  
  if ( ! db.FindXform( "image2.nii", "image1.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::cerr;
    return 1;
    }

  if ( xform != "xform12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::cerr;
    return 1;
    }
  
  if ( !inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::cerr;
    return 1;
    }

  return 0;
}

// test getting transformations when none actually exists
int
testImageXformDBFindXformNone()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  if ( ! db.AddImagePairXform( "xform12", true /*invertible*/, "image1.nii", "image2.nii" ) )
    return 1;
  
  std::string xform;
  bool inverse;
  
  if ( db.FindXform( "image1.nii", "image3.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup succeeded when it should have failed." << std::cerr;
    return 1;
    }

  if ( db.FindXform( "image3.nii", "image2.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup succeeded when it should have failed." << std::cerr;
    return 1;
    }

  if ( db.FindXform( "image3.nii", "image3.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup succeeded when it should have failed." << std::cerr;
    return 1;
    }
  
  return 0;
}

// test getting transformations when multiple different ones exist
int
testImageXformDBFindXformMultiple()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  if ( ! db.AddImagePairXform( "affine12", true /*invertible*/, "image1.nii", "image2.nii" ) || 
       ! db.AddRefinedXform( "nonrigid12", false /*invertible*/, "image1.nii", "image2.nii" ) )
    return 1;
  
  std::string xform;
  bool inverse;
  
  if ( ! db.FindXform( "image1.nii", "image2.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::cerr;
    return 1;
    }
  
  if ( xform != "nonrigid12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::cerr;
    return 1;
    }
  
  if ( inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::cerr;
    return 1;
    }

  if ( ! db.FindXform( "image2.nii", "image1.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::cerr;
    return 1;
    }
  
  if ( xform != "nonrigid12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::cerr;
    return 1;
    }
  
  if ( !inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::cerr;
    return 1;
    }

  return 0;
}

// test getting transformations when multiple refinement levels exist.
int
testImageXformDBFindXformLevels()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  if ( ! db.AddImagePairXform( "affine12", true /*invertible*/, "image1.nii", "image2.nii" ) || 
       ! db.AddRefinedXform( "nonrigid12", false /*invertible*/, "affine12" ) )
    return 1;
  
  std::string xform;
  bool inverse;
  
  if ( ! db.FindXform( "image1.nii", "image2.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::cerr;
    return 1;
    }
  
  if ( xform != "nonrigid12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::cerr;
    return 1;
    }
  
  if ( inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::cerr;
    return 1;
    }

  if ( ! db.FindXform( "image2.nii", "image1.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::cerr;
    return 1;
    }
  
  if ( xform != "nonrigid12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::cerr;
    return 1;
    }
  
  if ( !inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::cerr;
    return 1;
    }

  return 0;
}

// test getting transformations when multiple refinement levels exist and the refined transformation is based on an inverted initial transformation
int
testImageXformDBFindXformLevelsInverse()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  if ( ! db.AddImagePairXform( "affine12", true /*invertible*/, "image1.nii", "image2.nii" ) || 
       ! db.AddRefinedXform( "nonrigid21", false /*invertible*/, "affine12", true /*initInverse*/ ) )
    return 1;
  
  std::string xform;
  bool inverse;
  
  if ( ! db.FindXform( "image1.nii", "image2.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::cerr;
    return 1;
    }
  
  if ( xform != "affine12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::cerr;
    return 1;
    }
  
  if ( !inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::cerr;
    return 1;
    }

  if ( ! db.FindXform( "image2.nii", "image1.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::cerr;
    return 1;
    }
  
  if ( xform != "nonrigid21" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::cerr;
    return 1;
    }
  
  if ( inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::cerr;
    return 1;
    }

  return 0;
}

