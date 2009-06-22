/*
//
//  Copyright 1997-2004 Torsten Rohlfing
//  Copyright 2009 SRI International
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
#include <cmtkAffineRegistration.h>
#include <cmtkProtocolCallback.h>

#include <cmtkPerform.h>

/** Print some diagnostic information about the volume.
*/
void dumpVolume( const cmtkUniformVolume *volume )
{
  fprintf( stderr, "%d x %d x %d pixels, %f x %f x %f mm\n", 
	   volume->GetDims( CMTK_AXIS_X ), volume->GetDims( CMTK_AXIS_Y ),volume->GetDims( CMTK_AXIS_Z ), 
	   volume->Size[0], volume->Size[1], volume->Size[2] );
  
  const cmtk::TypedArray *data = volume->GetData();
  if ( data ) 
    {
    cmtk::Types::DataItem min, max;
    data->GetRange( min, max );
    fprintf( stderr, "data range: %f .. %f\n", (float) min, (float) max );
    }
}

const char *Investigator = "Rohlfing T";
const char *Site = "SRI International";
const char *Method = "Multi-resolution optimization of NMI";
const char *Date = __DATE__;
const char *PatientNumber = "undefined";
const char *FromModality = "undefined";
const char *ToModality = "undefined";

char setPatientNumber[128];
char setFromModality[128];
char setToModality[128];

void parseFilenames( const char* refFn, const char *fltFn )
{
#ifdef _WIN32
#define PATH_SEPARATOR '\\'
#else
#define PATH_SEPARATOR '/'
#endif
  const char *last, *previous;
  last = strrchr( refFn, PATH_SEPARATOR );
  if ( ! last ) return;
  for ( previous = last-1; *previous != PATH_SEPARATOR && previous != refFn; --previous );
  strncpy( setFromModality, previous+1, (last-previous-1) );
  last = previous;
  for ( previous = last-1; *previous != PATH_SEPARATOR && previous != refFn; --previous );
  strncpy( setPatientNumber, previous+1, (last-previous-1) );
  
  last = strrchr( fltFn, PATH_SEPARATOR );
  if ( ! last ) return;
  for ( previous = last-1; *previous != PATH_SEPARATOR && previous != fltFn; --previous );
  strncpy( setToModality, previous+1, (last-previous-1) );
  
  FromModality = setFromModality;
  ToModality = setToModality;
  PatientNumber = setPatientNumber;
}

/** Dump transformation in Vanderbilt format as transformed corner coordinates.
*/
void dumpTransformationVanderbilt( const cmtkAffineXform* affineXform,
				   const cmtkVector3D& size )
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
    v.XYZ[2] = k ? size.XYZ[2] : 0;
    for ( int j=0; j<2; ++j ) 
      {
      v.XYZ[1] = j ? size.XYZ[1] : 0;
      for ( int i=0; i<2; ++i, ++pointIdx ) 
	{
	v.XYZ[0] = i ? size.XYZ[0] : 0;
	cmtkVector3D w( v );
	affineXform->ApplyInPlace( w );
	fprintf( stdout, " %1d %10.4f %10.4f %10.4f %10.4f %11.4f %11.4f\n", pointIdx, v.XYZ[0], v.XYZ[1], v.XYZ[2], w.XYZ[0], w.XYZ[1], w.XYZ[2] );
	}
      }
    }
  fprintf( stdout, "\n\n-------------------------------------------------------------------------\n" );
}

void DoRegistration( const char* refFile, const char* fltFile )
{
  parseFilenames( refFile, fltFile );
  
// read first (reference) volume and dump diagnostic information.
  cmtk::UniformVolume::SmartPtr refVolume( cmtk::VolumeIO::Read( refFile ) );
  if ( !refVolume ) 
    {
    fprintf( stderr, "Could not read reference volume %s.\n", refFile );
    exit( 1 );
    }
  dumpVolume( refVolume );

// read second (floating) volume and dump diagnostic information.
  cmtk::UniformVolume::SmartPtr fltVolume( cmtk::VolumeIO::Read( fltFile ) );
  if ( !fltVolume ) 
    {
    fprintf( stderr, "Could not read floating volume %s.\n", fltFile );
    exit( 1 );
    }
  dumpVolume( fltVolume );
  
// open new scope so we can destruct automatic Registration instance and potentially check memory usage.
  {
// setup registration object with two volumes; as soon as volumes are
// linked to registration object, we can decrement their reference
// counters.
  cmtk::AffineRegistration Registration;
// Registration.SetMetric( 1 ); // MI
  Registration.SetVolume_1( refVolume );
  cmtkVector3D refSize( refVolume->Size );
  Registration.SetVolume_2( fltVolume );
  
// we want the centers of both images aligned intially.
  Registration.SetInitialAlignCenters( true );
  
// set optimization parameters to what they were during the original Vanderbilt submission.
  Registration.SetExploration( 8.0 );
  Registration.SetAccuracy( 0.01 );
  Registration.SetSampling( 1.0 );
  Registration.SetUseOriginalData( true );
  
// create a simple callback object that creates a detailed protocol of the
// registration process.
  cmtk::RegistrationCallback::SmartPtr callback( new cmtk::ProtocolCallback( "vanderbilt.txt" ) );
  Registration.SetCallback( callback );
  
// print estimated memory allocation
// fprintf( stderr, "Estimated peak memory allocation: %d Kbytes.\n",
// Registration.EstimatedMemoryAllocation() );
  
// run registration
  cmtk::CallbackResult result = Registration.Register();
  
// output diagnostic message regarding the final status of the
// registration.
  fprintf( stderr, "Registration " );
  switch ( result ) 
    {
    case CMTK_CALLBACK_OK :
      fprintf( stderr, "completed successfully.\n" ); break;
    case CMTK_CALLBACK_INTERRUPT :
      fprintf( stderr, "was interrupted by user.\n" ); break;
    case CMTK_CALLBACK_TIMEOUT :
      fprintf( stderr, "reached timeout.\n" ); break;
    case CMTK_CALLBACK_FAILED :
      fprintf( stderr, "failed miserably.\n" ); break;
    default:
      fprintf( stderr, "returned unknown status code.\n" ); break;
    }
  
// output transformation parameters.
  cmtk::AffineXform::SmartPtr affineXform = Registration.GetTransformation();
// fprintf( stderr, "\rResulting transformation parameters: \n" );
// for ( int idx = 0; idx < affineXform->ParamVectorDim(); ++idx )
// fprintf( stderr, "#%d: %f\n", idx,
// (float) affineXform->Parameters[idx] );
  
  dumpTransformationVanderbilt( affineXform->GetInverse(), refSize );
  
// Registration object is destructed here.
  }
}

int main( const int argc, const char *argv[] )
{
// say who we are
  const char *cp = strrchr ( argv[0], '/' );
  if ( ! cp )
    cp = argv[0];
  else
    ++cp;
  fprintf( stderr, "%s %s %s\n", cp, __DATE__, __TIME__ );

  double timeBaseline = cmtkPerform::GetTimeProcess();
  DoRegistration( argv[1], argv[2] );
  fprintf( stderr, "time: %f [s]\n", cmtk::Perform::GetTimeProcess() - timeBaseline );
  
  return 0;
}

