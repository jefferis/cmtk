/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2014 SRI International
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
  std::string inputPath;
  std::string outputPath;

  std::string outputLabelPath;
  std::string outputRigidPath;
  std::string outputAffinePath;

  cmtk::DetectPhantomMagphanEMR051::Parameters detectionParameters;
  
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Detect ADNI phantom landmarks in phantom image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool detects the locations of all spherical landmarks in a 3D image of the Magphan EMR051 structural imaging phantom (a.k.a. ADNI Phantom)." );

    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Detection", "Phantom Detection Options" );
    cl.AddSwitch( Key( "tolerant" ), &detectionParameters.m_TolerateTruncation, true, "Be tolerant of issues such as partially truncated marker spheres. "
		  "This should be used with caution only when necessary, and both the phantom image and detection results should be carefully inspected to identify the source of detection problems and verify reliable results." );
    cl.AddSwitch( Key( "any-orientation" ), &detectionParameters.m_StandardOrientation, false, "Do not assume standard orientation of the phantom, i.e., allow phantoms scanned upside-down. This makes detection of defective phantoms less robust." );
    cl.AddOption( Key( "erode-snr" ), &detectionParameters.m_ErodeSNR, "Erode SNR sphere by this distance prior to computing SNR estimate." );
    cl.AddOption( Key( "erode-cnr" ), &detectionParameters.m_ErodeCNR, "Erode each CNR sphere by this distance prior to computing CNR estimate." );
    cl.AddSwitch( Key( "refine-xform" ), &detectionParameters.m_RefineXformEachLandmark, true, "Refine estimated affine transformation after each new landmark is added." );
    cl.AddSwitch( Key( "refine-outliers" ), &detectionParameters.m_RefineOutliers, true, "Refine outlier landmarks based on estimated transformation after first sphere detection pass." );
    cl.AddSwitch( Key( "exclude-outliers" ), &detectionParameters.m_ExcludeOutliers, true, "Exclude outlier landmarks before fitting final transformations." );
    cl.AddSwitch( Key( "no-bias-correct-spheres" ), &detectionParameters.m_CorrectSphereBiasField, false, "Disable intensity bias field correction for each detected sphere. This will likely reduce accuracy of SNR/CNR estimates and also affect "
		  "localication accuracy of smaller spheres, but may be helpful in extreme cases where bias correction fails completely." )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( "write-labels" ), &outputLabelPath, "Output label image path. This image contains the mask of detected spheres, each labeled uniquely in their order in CMTK's ADNI phantom fiducial table." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddOption( Key( "write-rigid" ), &outputRigidPath, "Output path for the fitted rigid transformation from phantom space into RAS image standard space. This transformation defines where each sphere should be in the image." )
      ->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddOption( Key( "write-affine" ), &outputAffinePath, "Output path for the fitted affine transformation from phantom space into RAS image standard space. This is the closest linear-fit transformation, "
		  "and as such it includes scale and shear components not present in the fitted rigid transformations. Since these components are due to scanner miscalibration and distortion, this transformation DOES NOT "
		  "specify the correct spehre locations in the image, but rather, allows for quantification of scale miscalibration.")
      ->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
    cl.EndGroup();
    
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

  try
    {
    cmtk::DetectPhantomMagphanEMR051 detectionFilter( phantomImage, detectionParameters );

    // get expected and actual landmark locations
    cmtk::LandmarkList expectedLandmarks = detectionFilter.GetExpectedLandmarks();
    cmtk::LandmarkList detectedLandmarks = detectionFilter.GetDetectedLandmarks();

    try
      {
      // bring expected landmark locations from phantom image into physical space
      const cmtk::AffineXform phantomToPhysical( phantomImage->GetImageToPhysicalMatrix() );
      for ( cmtk::LandmarkList::Iterator it = expectedLandmarks.begin(); it != expectedLandmarks.end(); ++it )
	{
	it->m_Location = phantomToPhysical.Apply( it->m_Location );
	}
    
      // bring detected landmark locations from phantom image into physical space
      for ( cmtk::LandmarkList::Iterator it = detectedLandmarks.begin(); it != detectedLandmarks.end(); ++it )
	{
	it->m_Location = phantomToPhysical.Apply( it->m_Location );
	}
      }
    catch ( const cmtk::AffineXform::MatrixType::SingularMatrixException& )
      {
      cmtk::StdErr << "ERROR: singular image-to-physical space matrix; cannot map landmarks to image space.\n";
      throw cmtk::ExitException( 1 );
      }
      
    // match expected and detected landmarks
    cmtk::LandmarkPairList pairList( expectedLandmarks, detectedLandmarks );
    cmtk::DebugOutput( 2 ) << "INFO: detected and matched " << pairList.size() << " out of " << expectedLandmarks.size() << " expected landmarks.\n";
    
    if ( ! outputPath.empty() )
      cmtk::PhantomIO::Write( *(detectionFilter.GetDetectedPhantom()), outputPath );
    
    if ( ! outputLabelPath.empty() )
      cmtk::VolumeIO::Write( *(detectionFilter.GetDetectedSpheresLabelMap()), outputLabelPath );
    
    if ( ! outputAffinePath.empty() )
      cmtk::XformIO::Write( detectionFilter.GetPhantomToImageTransformationAffine(), outputAffinePath );
    
    if ( ! outputRigidPath.empty() )
      cmtk::XformIO::Write( detectionFilter.GetPhantomToImageTransformationRigid(), outputRigidPath );
    }
  catch ( const cmtk::DetectPhantomMagphanEMR051::OutsideFieldOfView& ex )
    {
    cmtk::StdErr << "ERROR: estimated location " << ex.m_Location << " puts landmark #" << ex.m_Idx << " (partly) outside image field of view.\n";
    throw cmtk::ExitException( 3 );
    }
  catch ( const cmtk::DetectPhantomMagphanEMR051::NoSphereInSearchRegion& )
    {
    cmtk::StdErr << "ERROR: unable to find sphere near expected location - most likely insufficient field of view.\n";
    throw cmtk::ExitException( 3 );
    }
  catch ( const cmtk::AffineXform::MatrixType::SingularMatrixException& )
    {
    cmtk::StdErr << "ERROR: singular matrix in cmtk::DetectPhantomMagphanEMR051 constructor\n";
    throw cmtk::ExitException( 1 );
    }

  return 0;
}

#include "cmtkSafeMain"
