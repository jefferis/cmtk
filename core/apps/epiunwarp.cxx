/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkSmartPtr.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformIO.h>

#include <Base/cmtkUnits.h>

#include <Registration/cmtkEchoPlanarUnwarpFunctional.h>

int
doMain
( const int argc, const char *argv[] )
{
  std::string inputImagePath1;
  std::string inputImagePath2;

  std::string outputImagePath1;
  std::string outputImagePath2;

  std::string outputDField;

  std::string writeJacobianPath1;
  std::string writeJacobianPath2;  

  byte phaseEncodeDirection = 1;
  bool flipPEPolar = true;
  bool initShiftCentersOfMass = true;

  double smoothnessConstraintWeight = 0;
  double foldingConstraintWeight = 0;
  int iterations = 10;

  double smoothSigmaMax = 8.0;
  double smoothSigmaMin = 0;
  double smoothSigmaDiff = 0.25;
  
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Unwarp Echo Planar Images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Correct B0 field inhomogeneity-induced distortion in Echo Planar Images (e.g., diffusion-weighted images) using two images acquired with opposing phase encoding directions." );

    typedef cmtk::CommandLine::Key Key;

    cl.BeginGroup( "Input", "Input Image Parameters" );    
    cmtk::CommandLine::EnumGroup<byte>::SmartPtr phaseEncodeGroup = cl.AddEnum( "phase-encode", &phaseEncodeDirection, "Define the phase-encoded image coordinate direction." );
    phaseEncodeGroup->AddSwitch( Key( 'y', "phase-encode-ap" ), 1, "Anterior/posterior phase encoding (this is the most common case)" );
    phaseEncodeGroup->AddSwitch( Key( 'z', "phase-encode-is" ), 2, "Top/bottom phase encoding (this is rare)" );
    phaseEncodeGroup->AddSwitch( Key( 'x', "phase-encode-lr" ), 0, "Lateral, left/right phase encoding (this is extremely rare)" );
    cl.AddSwitch( Key( "no-flip" ), &flipPEPolar, false, "Use this switch is the reverse phase-encoded image does not need to be flipped prior to unwarping. If normal and reverse phase-encoded image display in the same gross "
		  "orientation in 'triplanar', then flipping is not necessary and must be turned off using this switch." );
    cl.EndGroup();

    cl.BeginGroup( "Optimization", "Optimization Parameters" );    
    cl.AddSwitch( Key( "no-init-shift-com" ), &initShiftCentersOfMass, false, "Disable initialization of unwarping by shifting each row to align the centers of mass of forward and reverse acquisition. "
		  "Instead, use all-zero initial deformation field." );

    cl.AddOption( Key( "smooth-sigma-max" ), &smoothSigmaMax, "Maximum image smoothing kernel width for coarsest level of multi-scale computation." );
    cl.AddOption( Key( "smooth-sigma-min" ), &smoothSigmaMin, "Minimum image smoothing kernel width for finest level of multi-scale computation (0 = no smoothing; original image scale)." );
    cl.AddOption( Key( "smooth-sigma-diff" ), &smoothSigmaDiff, "Difference between image smoothing kernel widths between two successive levels of the multi-scale computation." );
    
    cl.AddOption( Key( "smoothness-constraint-weight" ), &smoothnessConstraintWeight, "Weight factor for the second-order smoothness constraint term in the unwarping cost function." );
    cl.AddOption( Key( "folding-constraint-weight" ), &foldingConstraintWeight, "Weight factor for the folding-prevention constraint term in the unwarping cost function." );
    cl.AddOption( Key( 'i', "iterations" ), &iterations, "Number of L-BFGS optimization iterations (per multi-scale level)." );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( "write-jacobian-fwd" ), &writeJacobianPath1, "Write Jacobian intensity correction map for forward image." );
    cl.AddOption( Key( "write-jacobian-rev" ), &writeJacobianPath2, "Write Jacobian intensity correction map for reverse-encoded image." );
    cl.EndGroup();

    cl.AddParameter( &inputImagePath1, "InputImage1", "First input image path - this is the standard b=0 image." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &inputImagePath2, "InputImage2", "Second input image path - this is the b=0 image with reversed phase encoding direction." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &outputImagePath1, "OutputImage1", "First output image path - this is the unwarped, corrected standard b=0 image." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddParameter( &outputImagePath2, "OutputImage2", "Second output image path - this is the unwarped, corrected reversed-encoding b=0 image." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddParameter( &outputDField, "OutputDField", "Path for deformation field (this can be applied to other images, e.g., diffusion-weighted images." )->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartPtr inputImage1 = cmtk::VolumeIO::ReadOriented( inputImagePath1 );
  cmtk::UniformVolume::SmartPtr inputImage2 = cmtk::VolumeIO::ReadOriented( inputImagePath2 );

  if ( flipPEPolar )
    inputImage2->ApplyMirrorPlane( phaseEncodeDirection );

  cmtk::EchoPlanarUnwarpFunctional func( inputImage1, inputImage2, phaseEncodeDirection, initShiftCentersOfMass );

  func.SetSmoothnessConstraintWeight( smoothnessConstraintWeight );
  func.SetFoldingConstraintWeight( foldingConstraintWeight );
  func.Optimize( iterations, cmtk::Units::GaussianSigma( smoothSigmaMax ), cmtk::Units::GaussianSigma( smoothSigmaMin ), cmtk::Units::GaussianSigma( smoothSigmaDiff ) );

  cmtk::VolumeIO::Write( *func.GetCorrectedImage( +1 ), outputImagePath1 );
  cmtk::VolumeIO::Write( *func.GetCorrectedImage( -1 ), outputImagePath2 );

  if ( !outputDField.empty() )
    {
    cmtk::DeformationField::SmartPtr dfield( func.GetDeformationField( +1 ) );
    cmtk::XformIO::Write( dfield, outputDField );
    }
    
  if ( !writeJacobianPath1.empty() )
    cmtk::VolumeIO::Write( *func.GetJacobianMap( +1 ), writeJacobianPath1 );

  if ( !writeJacobianPath2.empty() )
    cmtk::VolumeIO::Write( *func.GetJacobianMap( -1 ), writeJacobianPath2 );

  return 0;
}

#include "cmtkSafeMain"
