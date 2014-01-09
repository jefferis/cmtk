/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <System/cmtkCommandLine.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>

#include <IO/cmtkXformIO.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkWarpXform.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkMathUtil.h>

#include <stdio.h>
#include <math.h>

#include <limits>

cmtk::Types::Coordinate Accuracy = 0.01;
bool NoCheck = false;
bool Forward = false;

const char* StudyList = NULL;

// Greg's additions: option to calculate jacobian at given point
bool Jacobian = false;
bool NormalisedJacobian = false;
// return the scaling from the affine component only
bool ReturnGlobalScaling = false;
bool AffineOnly = false;
// read and write binary data
bool Binary = false;
// Input Files
const char* InputPath=NULL;
const char* OutputPath=NULL;

const char* FallbackInversePath = NULL;

FILE *outfile = stdout;
FILE *infile=stdin;

int
doMain( const int argc, const char *argv[] )
{
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Apply coordinate transformations to lists of point coordinates" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "THIS TOOL IS DEPRECATED. PLEASE USE streamxform INSTEAD.\n\n"
		       "This tool reads a list of 3D coordinates and applies a coordimnate transformation to them. The transformation can optionally be inverted. The transformed coordinates are then written to standard output." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "gregxform [options] transformation" );      

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( 'a', "accuracy" ), &Accuracy, "Approximation accuracy (maximum error)" );
    cl.AddSwitch( Key( 'x', "no-check" ), &NoCheck, true, "Disable accuracy checking and always output something" );
    cl.AddSwitch( Key( 'f', "forward" ), &Forward, true, "Apply forward transformation (default: backward)" );
    cl.AddOption( Key( 'F', "fallback-inverse" ), &FallbackInversePath, "Fallback inverse transformation (used for initialization and in case of failure)" );

    // Greg's additions:
    cl.AddOption( Key( 'i', "input-file" ), &InputPath, "Input path [default: STDIN]. ");
    cl.AddOption( Key( 'o', "output-file" ), &OutputPath, "Output path [default: STDOUT]. ");
	
    cl.AddSwitch( Key( 'n', "affine" ), &AffineOnly, true, "Apply affine transformation even if warp studylist specified" );
    cl.AddSwitch( Key( 'b', "binary" ), &Binary, true, "Read Binary input and produce binary output (floats)" );
    cl.AddSwitch( Key( 'j', "jacobian" ), &Jacobian, true, "Calculate jacobian determinant rather than mapping point (default: false)\
					\n\tNB: JDet is >1 if the sample volume is larger than the reference" );
    cl.AddSwitch( Key( 'J', "normalised-jacobian" ), &NormalisedJacobian, true, "Calculate jacobian determinant normalised by global scaling factor (default: false)" );
    cl.AddSwitch( Key( 'g', "return-global-scaling" ), &ReturnGlobalScaling, true, "Return global scaling factor ie the scaling due to initial affine" );

    cl.Parse( argc, argv );

    StudyList = cl.GetNext();
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }

  cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( StudyList ) );
  if ( ! xform )
    {
    cmtk::StdErr << "ERROR: could not read transformation\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::WarpXform::SmartPtr warpXform = cmtk::WarpXform::SmartPtr::DynamicCastFrom( xform );
  cmtk::AffineXform::SmartPtr affineXform = cmtk::AffineXform::SmartPtr::DynamicCastFrom( xform );
  cmtk::SplineWarpXform::SmartPtr splineWarp = cmtk::SplineWarpXform::SmartPtr ::DynamicCastFrom( warpXform );
  
  if ( !affineXform && warpXform )
    affineXform = warpXform->GetInitialAffineXform();
  
  if( !affineXform && AffineOnly )
    {
    cmtk::StdErr << "Unable to obtain affine transform from: "<<StudyList<<"\n" ;
    throw cmtk::ExitException( 2 );
    }
  
  cmtk::AffineXform::SmartPtr inverseAffineXform( NULL );
  try
    {
    if ( affineXform )
      inverseAffineXform = affineXform->GetInverse();
    }
  catch ( const cmtk::AffineXform::MatrixType::SingularMatrixException& )
    {
    cmtk::StdErr << "ERROR: singular matrix encountered in cmtk::AffineXform::GetInverse()\n";
    throw cmtk::ExitException( 1 );
    }

  const cmtk::Types::Coordinate globalScaling = (splineWarp) ? splineWarp->GetGlobalScaling() : affineXform->GetGlobalScaling();

  cmtk::Xform::SmartPtr fallbackInverseXform;
  if ( FallbackInversePath )
    {
    fallbackInverseXform = cmtk::Xform::SmartPtr( cmtk::XformIO::Read( FallbackInversePath ) );
    }
  
  // GJ: Now open files if required
  if(InputPath)
    {
    infile = fopen( InputPath, "r" );
    if(infile==NULL)  
      {
      cmtk::StdErr << "Unable to open input path: "<<InputPath<<"\n" ;
      throw cmtk::ExitException( 1 );
      }
    }
  if(OutputPath)
    {
    outfile = fopen( OutputPath, "w" );
    if(outfile==NULL)  
      {
      cmtk::StdErr << "Unable to open output path: "<<OutputPath<<"\n" ;
      throw cmtk::ExitException( 1 );
      }
    }	  
  
  if( ReturnGlobalScaling )
    {
    fprintf( outfile, "%lf\n", globalScaling );
    return 0;
    }
  
  while ( ! feof( infile ) ) 
    {
    float xyz[3];

    size_t numRead;
    if(Binary)
      {
      numRead = fread( xyz, sizeof(*xyz), 3, infile);
      } 
    else
      {
      char line[80];      
      fgets( line, 80, infile );
      if ( feof( infile ) )
	{
	break;
	}
      numRead = sscanf( line, "%20f %20f %20f", xyz, xyz+1, xyz+2 );
      }
    
    if ( numRead == 3 ) 
      {
      cmtk::FixedVector<3,cmtk::Types::Coordinate> v = cmtk::FixedVector<3,cmtk::Types::Coordinate>::FromPointer( xyz );
      cmtk::FixedVector<3,cmtk::Types::Coordinate> u( v );
      cmtk::FixedVector<3,cmtk::Types::Coordinate> uu;

      cmtk::Types::Coordinate error = 0;
      bool success = true;
      if ( splineWarp && ! AffineOnly ) 
	{
	if( Jacobian || NormalisedJacobian )
	  {
	  cmtk::Types::Coordinate j = splineWarp->GetJacobianDeterminant(u);
	  if (NormalisedJacobian)
	    j /= globalScaling;
	  
	  if(Binary)
	    {
	    fwrite( &j, sizeof(float), 1, outfile );
	    } 
	  else
	    { 
	    fprintf( outfile, "%f\n", j );
	    }
	  continue;
	  }
	
	if ( Forward )
	  {
	  v = splineWarp->Apply( v );
	  success = true;
	  }
	else
	  {
	  if ( fallbackInverseXform )
	    {
	    cmtk::Vector3D initialEstimate( fallbackInverseXform->Apply( v ) );
	    success = splineWarp->ApplyInverseWithInitial( v, v, initialEstimate, Accuracy );
	    }
	  else
	    {
	    success = splineWarp->ApplyInverse( v, v, Accuracy );
	    }
	  }
	if ( !success )
	  {
	  uu = splineWarp->Apply( v ) - u;
	  error = uu.RootSumOfSquares();

	  fprintf( stderr, "ERROR: %f %f %f is not inside target image or inversion failed (error = %f)\n", u[0], u[1], u[2], error );
	  }
	} 
      else
	{
	if ( Forward )
	  {
	  if ( affineXform )
	    v = affineXform->Apply( v );
	  }
	else
	  {
	  if ( inverseAffineXform )
	    v = inverseAffineXform->Apply( v );
	  }
	}
      if ( success || NoCheck ) 
	{
	const float outxyz[3]={(float) v[0],(float) v[1], (float) v[2]};
	if (Binary)
	  {
	  fwrite( outxyz, sizeof(float), 3, outfile );
	  if ( success )
	    {
	    cmtk::DebugOutput( 1 ).GetStream().printf( "%f %f %f\n",  outxyz[0], outxyz[1], outxyz[2] );
	    }
	  else
	    {
	    cmtk::DebugOutput( 1 ).GetStream().printf( "%f %f %f E %f\n",  outxyz[0], outxyz[1], outxyz[2], error );
	    }
	  } 
	else
	  {
	  if ( success )
	    {	   
	    fprintf( outfile, "%f %f %f\n", outxyz[0], outxyz[1], outxyz[2] );
	    }
	  else
	    {
	    fprintf( outfile, "%f %f %f E %f\n",  outxyz[0], outxyz[1], outxyz[2], error );
	    }
	  }
	} 
      else 
	{
	if(Binary)
	  {
	  const float nan3[3] = { std::numeric_limits<float>::signaling_NaN(), std::numeric_limits<float>::signaling_NaN(), std::numeric_limits<float>::signaling_NaN()};
	  fwrite( nan3, sizeof(float), 3, outfile );
	  }
	else
	  {
	  fputs( "ERR ERR ERR\n", outfile );
	  }
	}
      }
    }
  
  return 0;
}

#include "cmtkSafeMain"
