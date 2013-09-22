/*
//
//  Copyright 2010, 2013 SRI International
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

#include <Registration/cmtkImageXformDB.h>

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

  const cmtk::ImageXformDB::PrimaryKeyType key1 = db.FindImageSpaceID( "image1.nii" );
  if ( key1 == cmtk::ImageXformDB::NOTFOUND )
    {
    std::cerr << "No space key for image1." << std::endl;
    return 1;
    }

  db.AddImage( "image2.nii", "image1.nii" );

  const cmtk::ImageXformDB::PrimaryKeyType key2 = db.FindImageSpaceID( "image2.nii" );
  if ( key2 == cmtk::ImageXformDB::NOTFOUND )
    {
    std::cerr << "No space key for image2." << std::endl;
    return 1;
    }

  if ( key1 != key2 )
    {
    std::cerr << "Keys for image1 and image2 do not match." << std::endl;
    return 1;
    }

  db.AddImage( "image3.nii", "image2.nii" );

  const cmtk::ImageXformDB::PrimaryKeyType key3 = db.FindImageSpaceID( "image3.nii" );
  if ( key3 == cmtk::ImageXformDB::NOTFOUND )
    {
    std::cerr << "No space key for image3." << std::endl;
    return 1;
    }

  if ( key1 != key3 )
    {
    std::cerr << "Keys for image1 and image3 do not match." << std::endl;
    return 1;
    }


  db.AddImage( "image4.nii" );
  const cmtk::ImageXformDB::PrimaryKeyType key4 = db.FindImageSpaceID( "image4.nii" );
  if ( key4 == cmtk::ImageXformDB::NOTFOUND )
    {
    std::cerr << "No space key for image4." << std::endl;
    return 1;
    }

  if ( key4 == key1 )
    {
    std::cerr << "Keys for image1 and image4 match when they should not." << std::endl;
    return 1;
    }

  return 0;
}

// test adding an image to a database repeatedly without creating multiple entries.
int
testImageXformDBAddImageRepeat()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  db.AddImage( "image1.nii" );

  const cmtk::ImageXformDB::PrimaryKeyType key1 = db.FindImageSpaceID( "image1.nii" );
  if ( key1 == cmtk::ImageXformDB::NOTFOUND )
    {
    std::cerr << "No space key for image1 (first add)" << std::endl;
    return 1;
    }

  db.AddImage( "image1.nii" );

  const cmtk::ImageXformDB::PrimaryKeyType key2 = db.FindImageSpaceID( "image1.nii" );
  if ( key2 == cmtk::ImageXformDB::NOTFOUND )
    {
    std::cerr << "No space key for image1 (second add)" << std::endl;
    return 1;
    }

  if ( key1 != key2 )
    {
    std::cerr << "Two different keys." << std::endl;
    return 1;
    }

  const std::vector<std::string> list = db.GetSpaceImageList( db.FindImageSpaceID( "image1.nii" ) );
  if ( list.size() != 1 )
    {
    std::cerr << "Number of entries is not equal to 1." << std::endl;
    return 1;
    }

  return 0;
}

// test adding an image to a database
int
testImageXformDBAddImagePair()
{
  cmtk::ImageXformDB db( ":memory:" );
  db.DebugModeOn();
  db.AddImage( "image2.nii", "image1.nii" );

  const cmtk::ImageXformDB::PrimaryKeyType key1 = db.FindImageSpaceID( "image1.nii" );
  if ( key1 == cmtk::ImageXformDB::NOTFOUND )
    {
    std::cerr << "No space key for image1." << std::endl;
    return 1;
    }

  const cmtk::ImageXformDB::PrimaryKeyType key2 = db.FindImageSpaceID( "image2.nii" );
  if ( key2 == cmtk::ImageXformDB::NOTFOUND )
    {
    std::cerr << "No space key for image2." << std::endl;
    return 1;
    }

  if ( key1 != key2 )
    {
    std::cerr << "Keys for image1 and image2 do not match." << std::endl;
    return 1;
    }

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
    std::cerr << "DB add transformation." << std::endl;
    return 1;
    }

  if ( (db.FindXformLevel( "xform12" ) != 0) ||
       (db.FindXformLevel( "xform12refined" ) != 1) )
    {
    std::cerr << "DB transformation levels are incorrect." << std::endl;
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
    std::cerr << "DB transformation lookup failed." << std::endl;
    return 1;
    }

  if ( xform != "xform12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::endl;
    return 1;
    }

  if ( inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::endl;
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
    std::cerr << "DB transformation lookup failed." << std::endl;
    return 1;
    }

  if ( xform != "" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::endl;
    return 1;
    }
  
  if ( inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::endl;
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
    std::cerr << "DB transformation lookup failed." << std::endl;
    return 1;
    }

  if ( xform != "xform12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::endl;
    return 1;
    }
  
  if ( !inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::endl;
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
    std::cerr << "DB transformation lookup succeeded when it should have failed." << std::endl;
    return 1;
    }

  if ( db.FindXform( "image3.nii", "image2.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup succeeded when it should have failed." << std::endl;
    return 1;
    }

  if ( db.FindXform( "image3.nii", "image3.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup succeeded when it should have failed." << std::endl;
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
       ! db.AddImagePairXform( "nonrigid12", false /*invertible*/, "image1.nii", "image2.nii" ) )
    return 1;
  
  std::string xform;
  bool inverse;
  
  if ( ! db.FindXform( "image1.nii", "image2.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::endl;
    return 1;
    }
  
  if ( xform != "nonrigid12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::endl;
    return 1;
    }
  
  if ( inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::endl;
    return 1;
    }

  if ( ! db.FindXform( "image2.nii", "image1.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::endl;
    return 1;
    }
  
  if ( xform != "nonrigid12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::endl;
    return 1;
    }
  
  if ( !inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::endl;
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
    std::cerr << "DB transformation lookup failed." << std::endl;
    return 1;
    }
  
  if ( xform != "nonrigid12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::endl;
    return 1;
    }
  
  if ( inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::endl;
    return 1;
    }

  if ( ! db.FindXform( "image2.nii", "image1.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::endl;
    return 1;
    }
  
  if ( xform != "nonrigid12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::endl;
    return 1;
    }
  
  if ( !inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::endl;
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
    std::cerr << "DB transformation lookup failed." << std::endl;
    return 1;
    }
  
  if ( xform != "affine12" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::endl;
    return 1;
    }
  
  if ( inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::endl;
    return 1;
    }

  if ( ! db.FindXform( "image2.nii", "image1.nii", xform, inverse ) )
    {
    std::cerr << "DB transformation lookup failed." << std::endl;
    return 1;
    }
  
  if ( xform != "nonrigid21" )
    {
    std::cerr << "DB transformation returned wrong xform." << std::endl;
    return 1;
    }
  
  if ( inverse )
    {
    std::cerr << "DB transformation returned wrong inversion flag." << std::endl;
    return 1;
    }

  return 0;
}

