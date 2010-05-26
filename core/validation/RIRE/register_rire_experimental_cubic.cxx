 /*
//
//  Copyright 1997-2004 Torsten Rohlfing
//
//  Copyright 2009-2010 SRI International
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

#include <cmtkVolumeIO.h>
#include <cmtkImagePairAffineRegistration.h>
#include <cmtkAnatomicalOrientation.h>

#include <cmtkTimers.h>

#include <algorithm>

/** Print some diagnostic information about the volume.
*/
void dumpVolume( const cmtk::UniformVolume *volume )
{
  fprintf( stderr, "%d x %d x %d pixels, %f x %f x %f mm\n", 
	   volume->GetDims()[cmtk::AXIS_X], volume->GetDims()[cmtk::AXIS_Y], volume->GetDims()[cmtk::AXIS_Z], 
	   volume->Size[cmtk::AXIS_X], volume->Size[cmtk::AXIS_Y], volume->Size[cmtk::AXIS_Z] );
  
  const cmtk::TypedArray *data = volume->GetData();
  if ( data ) 
    {
    const cmtk::Types::DataItemRange range = data->GetRange();
    fprintf( stderr, "data range: %f .. %f\n", (float) range.m_LowerBound, (float) range.m_UpperBound );
    }
}

const char *Investigator = "Rohlfing T";
const char *Site = "SRI International";
const char *Method = "CMTK "CMTK_VERSION" -- multi-resolution optimization of NMI using cubic interpolation";
const char *Date = __DATE__;
const char *PatientNumber = "undefined";
const char *FromModality = "undefined";
const char *ToModality = "undefined";

char setPatientNumber[128];
char setFromModality[128];
char setToModality[128];

void
parseFilenames( const char* refFn, const char *fltFn )
{  
#ifdef WIN32
#define PATH_SEPARATOR '\\'
#else
#define PATH_SEPARATOR '/'
#endif

  const char *last, *previous;
  last = strrchr( refFn, PATH_SEPARATOR );
  if ( ! last ) return;

  for ( previous = last-1; *previous != PATH_SEPARATOR && previous != refFn; --previous ) {};
  strncpy( setFromModality, previous+1, std::min<int>( sizeof( setFromModality ) - 1, (last-previous-1) ) );

  last = previous;
  for ( previous = last-1; *previous != '_' && previous != refFn; --previous ) {};
  strncpy( setPatientNumber, previous+1, std::min<int>( sizeof( setPatientNumber ) - 1, (last-previous-1) ) );
  
  last = strrchr( fltFn, PATH_SEPARATOR );
  if ( ! last ) return;
  for ( previous = last-1; *previous != PATH_SEPARATOR && previous != fltFn; --previous ) {};
  strncpy( setToModality, previous+1, std::min<int>( sizeof( setToModality ) - 1, (last-previous-1) ) );
  
  FromModality = setFromModality;
  ToModality = setToModality;
  PatientNumber = setPatientNumber;
}

/** Dump transformation in Vanderbilt format as transformed corner coordinates.
*/
void 
dumpTransformationVanderbilt( const cmtk::AffineXform::MatrixType& matrix, const cmtk::Vector3D& size )
{
  fprintf( stdout, "-------------------------------------------------------------------------\n" );
  fprintf( stdout, "Transformation Parameters\n\n" );
  fprintf( stdout, "Investigator(s): %s\n", Investigator );
  fprintf( stdout, "Site: %s\n\n", Site );
  fprintf( stdout, "Method: %s\n", Method );
  fprintf( stdout, "Date: %s\n", Date);
  fprintf( stdout, "Patient number: %s\n", PatientNumber );
  fprintf( stdout, "From: %s\n", FromModality );
  fprintf( stdout, "To: %s\n\n", ToModality );
  
  fprintf( stdout, "Point x y z new_x new_y new_z\n\n" );
  
  cmtk::Vector3D v;
  int pointIdx = 1;
  for ( int k=0; k<2; ++k ) 
    {
    v[2] = k ? size[2] : 0;
    for ( int j=0; j<2; ++j ) 
      {
      v[1] = j ? size[1] : 0;
      for ( int i=0; i<2; ++i, ++pointIdx ) 
	{
	v[0] = i ? size[0] : 0;
	cmtk::Vector3D w( v );
	matrix.Multiply( w );
	fprintf( stdout, " %1d %10.4f %10.4f %10.4f %10.4f %11.4f %11.4f\n", pointIdx, v[0], v[1], v[2], w[0], w[1], w[2] );
	}
      }
    }
  fprintf( stdout, "\n\n-------------------------------------------------------------------------\n" );
}

void DoRegistration( const char* refFile, const char* fltFile )
{
  parseFilenames( refFile, fltFile );
  
// read first (reference) volume and dump diagnostic information.
  cmtk::UniformVolume::SmartPtr refVolume( cmtk::VolumeIO::ReadOriented( refFile ) );
  if ( !refVolume ) 
    {
    fprintf( stderr, "Could not read reference volume %s.\n", refFile );
    exit( 1 );
    }
  dumpVolume( refVolume );

// read second (floating) volume and dump diagnostic information.
  cmtk::UniformVolume::SmartPtr fltVolume( cmtk::VolumeIO::ReadOriented( fltFile ) );
  if ( !fltVolume ) 
    {
    fprintf( stderr, "Could not read floating volume %s.\n", fltFile );
    exit( 1 );
    }
  dumpVolume( fltVolume );

// open new scope so we can destruct automatic Registration instance and potentially check memory usage.
  {
// setup registration object with two volumes
  cmtk::ImagePairAffineRegistration Registration;
  Registration.SetVolume_1( refVolume );
  cmtk::Vector3D refSize( refVolume->Size );
  Registration.SetVolume_2( fltVolume );
  Registration.SetFloatingImageInterpolation( cmtk::Interpolators::CUBIC );
  
// we want the centers of both images aligned intially.
// we want the centers of mass of both images aligned intially.
  Registration.SetInitializer( cmtk::MakeInitialAffineTransformation::COM );
  
// set optimization parameters to what they were during the original Vanderbilt submission.
  Registration.m_MaxStepSize = 8.0;
  Registration.m_MinStepSize = 0.01;
  Registration.m_Sampling = 1.0;
  Registration.m_UseOriginalData = true;

// run registration
  cmtk::CallbackResult result = Registration.Register();
  
// output diagnostic message regarding the final status of the registration.
  fprintf( stderr, "Registration " );
  switch ( result ) 
    {
    case cmtk::CALLBACK_OK :
      fprintf( stderr, "completed successfully.\n" ); break;
    case cmtk::CALLBACK_INTERRUPT :
      fprintf( stderr, "was interrupted by user.\n" ); break;
    case cmtk::CALLBACK_TIMEOUT :
      fprintf( stderr, "reached timeout.\n" ); break;
    case cmtk::CALLBACK_FAILED :
      fprintf( stderr, "failed miserably.\n" ); break;
    default:
      fprintf( stderr, "returned unknown status code.\n" ); break;
    }
  
// output transformation.
  const cmtk::AffineXform::MatrixType refMatrix = refVolume->GetImageToPhysicalMatrix ();
  const cmtk::AffineXform::MatrixType fltMatrix = fltVolume->GetImageToPhysicalMatrix ();
  cmtk::AffineXform::SmartPtr affineXform = Registration.GetTransformation();

  cmtk::AffineXform::MatrixType concatMatrix = refMatrix;
  (concatMatrix.Invert() *= affineXform->Matrix) *= fltMatrix;
  dumpTransformationVanderbilt( concatMatrix, refSize );
  
// Registration object is destructed here.
  }
}

int
main
( const int argc, const char *argv[] )
{
// say who we are
  const char *cp = strrchr ( argv[0], '/' );
  if ( ! cp )
    cp = argv[0];
  else
    ++cp;
  fprintf( stderr, "%s %s %s\n", cp, __DATE__, __TIME__ );

  if ( argc != 3 )
    {
    fputs( "ERROR: this tool needs exactly two parameters, the reference and floating image paths\n", stderr );
    exit( 1 );
    }

  double timeBaseline = cmtk::Timers::GetTimeProcess();
  DoRegistration( argv[1], argv[2] );
  fprintf( stderr, "time: %f [s]\n", cmtk::Timers::GetTimeProcess() - timeBaseline );
  
  return 0;
}

