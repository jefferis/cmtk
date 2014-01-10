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
//  $Revision: 5129 $
//
//  $LastChangedDate: 2014-01-09 13:48:15 -0800 (Thu, 09 Jan 2014) $
//
//  $LastChangedBy: torstenrohlfing $
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

#include <Base/cmtkDetectedPhantomMagphanEMR051.h>
#include <Base/cmtkLandmarkList.h>
#include <Base/cmtkLandmarkPairList.h>
#include <Base/cmtkFitPolynomialToLandmarks.h>

#include <vector>

int
doMain( const int argc, const char* argv[] )
{
  std::string inputPhantomPath;
  std::string inputImagePath;

  cmtk::Types::Coordinate residualThreshold = 5.0;

  bool fitInverse = false;
  byte degree = 4;

  std::string outputXform;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Unwarp T1-weighted MR image using a polynomial transformation based on a phantom description" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool computes a  polynomial transformation to unwarp an image. The transformation is based on expected and detected landmarks in an image of a structural phantom "
		       "acquired on the same scanner. Use the 'detect_adni_phantom' tool to detect landmarks of the ADNI Phantom in an image and generate a phantom description file suitable for use with this tool." );

    typedef cmtk::CommandLine::Key Key;    
    cl.BeginGroup( "Fitting", "Fitting Options" );
    cl.AddOption( Key( "degree" ), &degree, "Degree of the fitted polynomial transformation." );
    cl.AddSwitch( Key( "fit-inverse" ), &fitInverse, true, "Fit inverse transformation - this is useful for computing a Jacobian volume correction map (using 'reformatx') without having to numerically invert the fitted unwarping "
		  "transformation." );
    cl.EndGroup();

    cl.AddParameter( &inputPhantomPath, "InputPhantom", "Input path of the XML file describing a phantom previously detected in an image." );
    cl.AddParameter( &inputImagePath, "InputImage", "Input image path. This is the image that is unwarped. It is important that this image be acquired on the same scanner (not only the same model but the very machine) "
		     "on which the phantom image was also acquired, preferably in close temporal proximity. Also, both this and the phantom image must share and specify the same physical image coordinates, i.e., only images in "
		     "NIFTI or NRRD format can be used." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &outputXform, "OutputXform", "Output transformation path. This is the affine phantom-to-image coordinate transformation fitted to the detected landmark spheres." )
      ->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  // read phantom description
  cmtk::DetectedPhantomMagphanEMR051::SmartPtr phantom( cmtk::PhantomIO::Read( inputPhantomPath ) );
  cmtk::DebugOutput( 5 ) << "INFO: read phantom with " << phantom->LandmarkPairsList().size() << " landmarks.\n";  

  cmtk::UniformVolume::SmartConstPtr unwarpImage = cmtk::VolumeIO::ReadOriented( inputImagePath );
  try
    {
    phantom->ApplyXformToLandmarks( cmtk::AffineXform( unwarpImage->GetImageToPhysicalMatrix().GetInverse() ) );
    }
  catch ( const cmtk::AffineXform::MatrixType::SingularMatrixException& )
    {
    cmtk::StdErr << "ERROR: singular image-to-physical space matrix cannot be inverted.\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::LandmarkPairList pairList;
  for ( std::list<cmtk::LandmarkPair>::const_iterator it = phantom->LandmarkPairsList().begin(); it != phantom->LandmarkPairsList().end(); ++it )
    {
    if ( it->m_Precise ) // exclude all unprecise landmarks
      {
      if ( it->m_Residual < residualThreshold ) // exclude outliers based on residual
	{
	if ( fitInverse )
	  {
	  pairList.push_back( it->GetSwapSourceTarget() );
	  }
	else
	  {
	  pairList.push_back( *it );
	  }
	}
      }
    }

  cmtk::DebugOutput( 2 ) << "INFO: using " << pairList.size() << " out of " << phantom->LandmarkPairsList().size() << " total phantom landmarks as fiducials.\n";
  
  // fit polynomial transformation to landmark pairs.
  cmtk::PolynomialXform::SmartConstPtr polyXform = cmtk::FitPolynomialToLandmarks( pairList, degree ).GetPolynomialXform();
  
  // writing resulting transformation
  cmtk::XformIO::Write( polyXform, outputXform );

  return 0;
}

#include "cmtkSafeMain"
