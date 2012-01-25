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
//  $Revision: 326 $
//
//  $LastChangedDate: 2009-07-29 15:08:42 -0700 (Wed, 29 Jul 2009) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>
#include <System/cmtkProgressConsole.h>

#include <IO/cmtkVolumeIO.h>
#include <Segmentation/cmtkAtlasSegmentation.h>

#include <list>
#include <string>

int
doMain( const int argc, const char* argv[] )
{
  bool fast = false;
  
  const char* targetImageName = NULL;
  const char* outImageName = NULL;

  std::string channelSRI24( "spgr" );
  std::string labelsSRI24( "tzo116plus" );
  
  try
    {
    cmtk::CommandLine  cl( cmtk::CommandLine::PROPS_XML );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Segmentation using SRI24 atlas" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Register a target image to a selected channel of the SRI24 human brain atlas, then reformat one of the atlas label maps to the target image. "
		       "Note: it is assume that the target image is skull-stripped, i.e., contains only the brain." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Segmentation" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'f', "fast" ), &fast, true, "Fast mode." );

    cmtk::CommandLine::EnumGroup<std::string>::SmartPtr channelGroup = 
      cl.AddEnum( "registration-channel", &channelSRI24, "The SRI24 channel used for registration to the target image." );
    channelGroup->AddSwitch( Key( "spgr" ), "spgr", "SPGR (T1-weighted) structural channel" );
    channelGroup->AddSwitch( Key( "early-fse" ), "erly", "Early-echo (PD-weighted) fast spin echo channel" );
    channelGroup->AddSwitch( Key( "late-fse" ), "late", "Late-echo (T2-weighted) fast spin echo channel" );
    channelGroup->AddSwitch( Key( "fa" ), "fa", "Fractional anisotropy channel, derived from diffusion tensor images" );
    
    cmtk::CommandLine::EnumGroup<std::string>::SmartPtr labelsGroup = 
      cl.AddEnum( "label-map", &labelsSRI24, "The SRI24 label map that is reformatted to the target image." );
    labelsGroup->AddSwitch( Key( "tzo116plus" ), "tzo116plus", "Extended cortical parcellation template based on Tzourio-Mazoyer AAL 116 region template." );
    labelsGroup->AddSwitch( Key( "lpba40" ), "lpba40", "Template based on the 40 subject LONI Probabilistic Brain Atlas segmentation." );
    labelsGroup->AddSwitch( Key( "tissue" ), "tissue", "SRI24 maximum likelihood three compartment (GM, WM, CSF) tissue segmentation map." );
    
    cl.AddParameter( &targetImageName, "TargetImage", "Target image path. This is the image to be segmented." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &outImageName, "OutputImage", "Output image path. This is where the reformatted label map is written." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_LABELS | cmtk::CommandLine::PROPS_OUTPUT );

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }
  
  // Instantiate programm progress indicator.
  cmtk::ProgressConsole progressIndicator( "Segmentation Using the SRI24 Atlas" );

  cmtk::UniformVolume::SmartPtr targetImg( cmtk::VolumeIO::ReadOriented( targetImageName ) );
  if ( !targetImg ) 
    {
    cmtk::StdErr << "ERROR: could not read target image " << targetImageName << "\n";
    throw cmtk::ExitException( 1 );
    }
  
  std::string atlasImageName = std::string( CMTK_ROOT_PATH_SRI24 ) + "/" + channelSRI24 + ".nii";
  cmtk::UniformVolume::SmartPtr atlasImg( cmtk::VolumeIO::ReadOriented( atlasImageName.c_str() ) );
  if ( !atlasImg ) 
    {
    cmtk::StdErr << "ERROR: could not read atlas image " << atlasImageName << "\n";
    throw cmtk::ExitException( 1 );
    }
  
  std::string atlasLabelName = std::string( CMTK_ROOT_PATH_SRI24 ) + "/" + labelsSRI24 + ".nii";
  cmtk::UniformVolume::SmartPtr atlasLbl( cmtk::VolumeIO::ReadOriented( atlasLabelName.c_str() ) );
  if ( !atlasLbl ) 
    {
    cmtk::StdErr << "ERROR: could not read atlas labels " << atlasLabelName << "\n";
    throw cmtk::ExitException( 1 );
    }
    
  cmtk::AtlasSegmentation segment( targetImg, atlasImg, atlasLbl );
  segment.SetFast( fast );

  cmtk::VolumeIO::Write( *(segment.GetLabelMap()), outImageName );

  return 0;
}

#include "cmtkSafeMain"
