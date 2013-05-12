/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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
#include <System/cmtkExitException.h>
#include <System/cmtkDebugOutput.h>

#include <IO/cmtkPhantomIO.h>
#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformIO.h>

#include <Segmentation/cmtkDetectPhantomMagphanEMR051.h>

#include <vector>

int
doMain( const int argc, const char* argv[] )
{
  const char* inputPath = NULL;
  const char* outputPath = NULL;

  const char* outputLabelPath = NULL;
  const char* outputRigidPath = NULL;
  const char* outputAffinePath = NULL;

  bool tolerant = false;
  byte erodePixelsSNR = 2;
  byte erodePixelsCNR = 2;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Detect ADNI phantom landmarks in phantom image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool detects the locations of all spherical landmarks in a 3D image of the Magphan EMR051 structural imaging phantom (a.k.a. ADNI Phantom)." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( "tolerant" ), &tolerant, true, "Be tolerant of issues such as partially truncated marker spheres. This should be considered a last-ditch resort, and both phantom image and detection results should be carefully inspected." )
      cl.AddOption( Key( "erode-pixels-snr" ), &erodePixelsSNR, "Erode SNR sphere by this many pixels prior to computing SNR estimate." );
      cl.AddOption( Key( "erode-pixels-cnr" ), &erodePixelsSNR, "Erode each CNR sphere by this many pixels prior to computing CNR estimate." );
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddOption( Key( "write-labels" ), &outputLabelPath, "Output label image path. This image contains the mask of detected spheres, each labeled uniquely in their order in CMTK's ADNI phantom fiducial table." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddOption( Key( "write-rigid" ), &outputRigidPath, "Output path for the fitted rigid transformation from phantom space into RAS image standard space. This transformation defines where each sphere should be in the image." )
      ->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddOption( Key( "write-affine" ), &outputAffinePath, "Output path for the fitted affine transformation from phantom space into RAS image standard space. This is the closest linear-fit transformation, "
		  "and as such it includes scale and shear components not present in the fitted rigid transformations. Since these components are due to scanner miscalibration and distortion, this transformation DOES NOT "
		  "specify the correct spehre locations in the image, but rather, allows for quantification of scale miscalibration.")
      ->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.AddParameter( &inputPath, "InputImage", "Input image path. This is the image in which spheres are detected." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &outputPath, "OutputXML", "Output path for the XML file describing the dected phantom." )->SetProperties( cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  cmtk::UniformVolume::SmartPtr phantomImage( cmtk::VolumeIO::ReadOriented( inputPath ) );
  const cmtk::AffineXform phantomToPhysical( phantomImage->GetImageToPhysicalMatrix() );

  try
    {
    cmtk::DetectPhantomMagphanEMR051 detectionFilter( phantomImage, tolerant );

    // get expected landmark locations
    cmtk::LandmarkList expectedLandmarks = detectionFilter.GetExpectedLandmarks();
    // bring expected landmark locations from phantom image into physical space
    for ( cmtk::LandmarkList::Iterator it = expectedLandmarks.begin(); it != expectedLandmarks.end(); ++it )
      {
      it->m_Location = phantomToPhysical.Apply( it->m_Location );
      }
    
    // get detected landmark locations
    cmtk::LandmarkList actualLandmarks = detectionFilter.GetDetectedLandmarks();
    // bring detected landmark locations from phantom image into physical space
    for ( cmtk::LandmarkList::Iterator it = expectedLandmarks.begin(); it != expectedLandmarks.end(); ++it )
      {
      it->m_Location = phantomToPhysical.Apply( it->m_Location );
      }
    
    // match expected and detected landmarks
    cmtk::LandmarkPairList pairList( expectedLandmarks, actualLandmarks );
    cmtk::DebugOutput( 2 ) << "INFO: detected and matched " << pairList.size() << " out of " << expectedLandmarks.size() << " expected landmarks.\n";
    
    if ( outputPath )
      cmtk::PhantomIO::Write( *(detectionFilter.GetDetectedPhantom()), outputPath );
    
    if ( outputLabelPath )
      cmtk::VolumeIO::Write( *(detectionFilter.GetDetectedSpheresLabelMap()), outputLabelPath );
    
    if ( outputAffinePath )
      cmtk::XformIO::Write( detectionFilter.GetPhantomToImageTransformationAffine(), outputAffinePath );
    
    if ( outputRigidPath )
      cmtk::XformIO::Write( detectionFilter.GetPhantomToImageTransformationRigid(), outputRigidPath );
    }
  catch ( const cmtk::DetectPhantomMagphanEMR051::OutsideFieldOfView& ex )
    {
    cmtk::StdErr << "ERROR: estimated location " << ex.m_Location << " puts landmark #" << ex.m_Idx << " (partly) outside image field of view.\n";
    return 3;
    }
  catch ( const cmtk::DetectPhantomMagphanEMR051::NoSphereInSearchRegion )
    {
    cmtk::StdErr << "ERROR: unable to find sphere near expected location - most likely insufficient field of view.\n";
    return 3;
    }

  return 0;
}

#include "cmtkSafeMain"
