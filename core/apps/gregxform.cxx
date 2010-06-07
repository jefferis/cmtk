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

#include <cmtkCommandLine.h>
#include <cmtkConsole.h>

#include <cmtkXformIO.h>

#include <cmtkAffineXform.h>
#include <cmtkWarpXform.h>
#include <cmtkSplineWarpXform.h>
#include <cmtkMathUtil.h>

#include <cstdio>
#include <cmath>

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace gregxform
{
#endif
cmtk::Types::Coordinate Accuracy = 0.01;
bool NoCheck = false;
bool Forward = false;

const char* StudyList = NULL;

bool Verbose = false;

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
main( int argc, char *argv[] )
{
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Apply coordinate transformations to lists of point coordinates" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] transformation" );      

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
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Print each point to STDERR (as well as stdout)" );

    cl.Parse();

    StudyList = cl.GetNext();
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    exit( 1 );
    }

  cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( StudyList, Verbose ) );
  if ( ! xform )
    {
    cmtk::StdErr << "ERROR: could not read transformation\n";
    exit( 1 );
    }

  cmtk::WarpXform::SmartPtr warpXform = cmtk::WarpXform::SmartPtr::DynamicCastFrom( xform );
  cmtk::AffineXform::SmartPtr affineXform = cmtk::AffineXform::SmartPtr::DynamicCastFrom( xform );
  cmtk::SplineWarpXform::SmartPtr splineWarp = cmtk::SplineWarpXform::SmartPtr ::DynamicCastFrom( warpXform );
  
  if ( !affineXform && warpXform )
    affineXform = warpXform->GetInitialAffineXform();
  
  if( !affineXform && AffineOnly )
    {
    cmtk::StdErr << "Unable to obtain affine transform from: "<<StudyList<<"\n" ;
    exit( 2 );
    }
  
  cmtk::AffineXform::SmartPtr inverseAffineXform( NULL );
  if ( affineXform )
    inverseAffineXform = affineXform->GetInverse();

  const float globalScaling = (splineWarp) ? splineWarp->GetGlobalScaling() : affineXform->GetGlobalScaling();

  cmtk::Xform::SmartPtr fallbackInverseXform;
  if ( FallbackInversePath )
    {
    fallbackInverseXform = cmtk::Xform::SmartPtr( cmtk::XformIO::Read( FallbackInversePath, Verbose ) );
    }
  
  // GJ: Now open files if required
  if(InputPath)
    {
    infile = fopen( InputPath, "r" );
    if(infile==NULL)  
      {
      cmtk::StdErr << "Unable to open input path: "<<InputPath<<"\n" ;
      exit( 1 );
      }
    }
  if(OutputPath)
    {
    outfile = fopen( OutputPath, "w" );
    if(outfile==NULL)  
      {
      cmtk::StdErr << "Unable to open output path: "<<OutputPath<<"\n" ;
      exit( 1 );
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
      numRead = sscanf( line, "%f %f %f", xyz, xyz+1, xyz+2 );
      }
    
    if ( numRead == 3 ) 
      {
      cmtk::FixedVector<3,cmtk::Types::Coordinate> v( xyz );
      cmtk::FixedVector<3,cmtk::Types::Coordinate> u( v );
      cmtk::FixedVector<3,cmtk::Types::Coordinate> uu;

      float error = 0;
      bool success = true;
      if ( splineWarp && ! AffineOnly ) 
	{
	if( Jacobian || NormalisedJacobian )
	  {
	  float j = splineWarp->GetJacobianDeterminant(u);
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
	  splineWarp->ApplyInPlace( v );
	  success = true;
	  }
	else
	  {
	  if ( fallbackInverseXform )
	    {
	    cmtk::Vector3D initialEstimate( fallbackInverseXform->Apply( v ) );
	    success = splineWarp->ApplyInverseInPlaceWithInitial( v, initialEstimate, Accuracy );
	    }
	  else
	    {
	    success = splineWarp->ApplyInverseInPlace( v, Accuracy );
	    }
	  }
	if ( !success )
	  {
	  uu = v;
	  splineWarp->ApplyInPlace( uu );
	  uu -= u;
	  error = uu.RootSumOfSquares();

	  fprintf( stderr, "ERROR: %f %f %f is not inside target image or inversion failed (error = %f)\n", u[0], u[1], u[2], error );
	  }
	} 
      else
	{
	if ( Forward )
	  {
	  if ( affineXform )
	    affineXform->ApplyInPlace( v );
	  }
	else
	  {
	  if ( inverseAffineXform )
	    inverseAffineXform->ApplyInPlace( v );
	  }
	}
      if ( success || NoCheck ) 
	{
	const float outxyz[3]={(float) v[0],(float) v[1], (float) v[2]};
	if (Verbose) 
	  {
	  if ( success )
	    {
	    fprintf( stderr, "%f %f %f\n",  outxyz[0], outxyz[1], outxyz[2] );
	    }
	  else
	    {
	    fprintf( stderr, "%f %f %f E %f\n",  outxyz[0], outxyz[1], outxyz[2], error );
	    }
	  }
	if (Binary)
	  {
	  fwrite( outxyz, sizeof(float), 3, outfile );
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
	  const float nan3[3] = { cmtk::MathUtil::GetFloatNaN(), cmtk::MathUtil::GetFloatNaN(), cmtk::MathUtil::GetFloatNaN()};
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
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace gregxform
} // namespace apps
} // namespace cmtk
#endif

