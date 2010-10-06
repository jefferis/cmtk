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

#include <System/cmtkConsole.h>
#include <System/cmtkCommandLine.h>

#include <Base/cmtkMetaInformationObject.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkTransformChangeToSpaceAffine.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformIO.h>
#include <IO/cmtkAffineXformITKIO.h>

int 
main( const int argc, const char* argv[] )
{
  bool verbose = false;
  
  const char* inputPath = NULL;
  const char* outputPath = NULL;

  const char* fixedImagePath = NULL;
  const char* movingImagePath = NULL;
  
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Convert affine transformations to ITK format and correct for differences in image coordinate conventions" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &verbose, true, "Verbose mode" );

    cl.BeginGroup( "Output", "Output parameters" )->SetProperties( cmtk::CommandLine::PROPS_NOXML );
    cl.AddOption( Key( "fixed-image" ), &fixedImagePath, "Override transformation's fixed/reference image (if any) with this one." );
    cl.AddOption( Key( "moving-image" ), &movingImagePath, "Override transformation's moving/floating image (if any) with this one." );
    cl.EndGroup();
    
    cl.AddParameter( &inputPath, "InputPath", "CMTK input transformation path" )->SetProperties( cmtk::CommandLine::PROPS_XFORM );
    cl.AddParameter( &outputPath, "OutputPath", "ITK output transformation path" )->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.Parse( argc, argv );

    }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  cmtk::AffineXform::SmartConstPtr xform = cmtk::AffineXform::SmartConstPtr::DynamicCastFrom( cmtk::XformIO::Read( inputPath, verbose ) );

  if ( !fixedImagePath )
    {
    fixedImagePath = xform->GetMetaInfo( cmtk::META_XFORM_FIXED_IMAGE_PATH ).c_str();
    }

  if ( !movingImagePath )
    {
    movingImagePath = xform->GetMetaInfo( cmtk::META_XFORM_MOVING_IMAGE_PATH ).c_str();
    }

  cmtk::UniformVolume::SmartConstPtr fixedImage = cmtk::VolumeIO::ReadGridOriented( fixedImagePath, verbose );
  if ( ! fixedImage )
    {
    cmtk::StdErr << "ERROR: could not read fixed image " << fixedImagePath << "\n";
    exit( 1 );
    }

  cmtk::UniformVolume::SmartConstPtr movingImage = cmtk::VolumeIO::ReadGridOriented( movingImagePath, verbose );
  if ( ! movingImage )
    {
    cmtk::StdErr << "ERROR: could not read moving image " << movingImagePath << "\n";
    exit( 1 );
    }

  cmtk::TransformChangeToSpaceAffine toNative( *(xform), *(fixedImage), *(movingImage) );
  cmtk::AffineXformITKIO::Write( outputPath, toNative.GetTransformation() );
  
  return 0;
}

