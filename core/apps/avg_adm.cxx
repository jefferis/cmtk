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

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>
#include <System/cmtkProgressConsole.h>
#include <System/cmtkThreads.h>
#include <System/cmtkTimers.h>
#include <System/cmtkStrUtility.h>

#include <IO/cmtkClassStream.h>
#include <IO/cmtkTypedStreamStudylist.h>
#include <IO/cmtkVolumeIO.h>

#include <Base/cmtkActiveDeformationModel.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkFilterMask.h>

#include <Registration/cmtkReformatVolume.h>

#include <string.h>
#include <map>
#include <string>
#include <vector>
#include <list>
#include <iostream>

#ifdef CMTK_BUILD_MPI
#  include <mpi.h>
#endif

bool Verbose = false;

bool Labels = false;
bool Jacobian = false;

bool AutoScale = false;

cmtk::ScalarDataType UserDataType = cmtk::TYPE_NONE;
cmtk::Interpolators::InterpolationEnum Interpolation = cmtk::Interpolators::LINEAR;

const char* ReplaceFrom;
std::map<std::string,std::string> ReplaceMap;

void AddReplaceFrom( const char* arg )
{
  ReplaceFrom = arg;
}

void AddReplaceTo( const char* arg )
{
  ReplaceMap[std::string(ReplaceFrom)] = std::string( arg );
}

const char* OutImageName = NULL;
const char* OutWarpName = NULL;
const char* InWarpName = NULL;

bool IncludeScale = true;
bool IncludeReferenceModel = true;
bool IncludeReferenceData = true;

bool WriteIncludeAffine = false;

float PaddingValue = 0;
bool SetPaddingValue = false;

float NewPaddingValue = 0;
bool ReplacePaddingValue = false;

cmtk::Types::DataItem PadOutValue = 0;
bool HavePadOutValue = false;

int
doMain( int argc, char** argv )
{
#ifdef CMTK_BUILD_MPI
  MPI::Init( argc, argv );
  const int mpiRank = MPI::COMM_WORLD.Get_rank();
  const int mpiSize = MPI::COMM_WORLD.Get_size();
#endif

  std::list<const char*> inFileList;

  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Average using ADM" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Compute average-shape average-intensity images and deformation maps using an active deformation model." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] list0 [list1 ...]" );
    
    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddSwitch( Key( 'l', "label" ), &Labels, true, "Label mode (as opposed to grey)" );
    cl.AddOption( Key( "set-padding" ), &PaddingValue, "Set padding value in input files", &SetPaddingValue );
    cl.AddOption( Key( "replace-padding" ), &NewPaddingValue, "Replace padding value in input files", &ReplacePaddingValue );
    cl.AddOption( Key( 'p', "pad-out" ), &PadOutValue, "Padding value for output file", &HavePadOutValue );
    cl.AddSwitch( Key( 'j', "jacobian" ), &Jacobian, true, "Average Jacobian determinants of deformations." );
    cl.AddSwitch( Key( 'a', "auto-scale" ), &AutoScale, true, "Auto-scale image intensities" );

    cl.AddSwitch( Key( "char" ), &UserDataType, cmtk::TYPE_CHAR, "Output is 8 bits, signed [default: automatic]" );
    cl.AddSwitch( Key( "byte" ), &UserDataType, cmtk::TYPE_BYTE, "Output is 8 bits, unsigned" );
    cl.AddSwitch( Key( "short" ), &UserDataType, cmtk::TYPE_SHORT, "Output is 16 bits, signed" );
    cl.AddSwitch( Key( "ushort" ), &UserDataType, cmtk::TYPE_USHORT, "Output is 16 bits, unsigned" );
    cl.AddSwitch( Key( "int" ), &UserDataType, cmtk::TYPE_INT, "Output is 32 bits signed" );
    cl.AddSwitch( Key( "float" ), &UserDataType, cmtk::TYPE_FLOAT, "Output is 32 bits floating point" );
    cl.AddSwitch( Key( "double" ), &UserDataType, cmtk::TYPE_DOUBLE, "Output is 64 bits floating point\n" );

    cl.AddSwitch( Key( 'S', "no-scale-model" ), &IncludeScale, false, "Exclude scale from statistical deformation model" );
    cl.AddSwitch( Key( 'R', "no-ref-model" ), &IncludeReferenceModel, false, "Exclude reference coordinate frame from deformation model" );
    cl.AddSwitch( Key( 'r', "no-ref-data" ), &IncludeReferenceData, false, "Exclude reference data from averaging" );

    cl.AddOption( Key( 'o', "output" ), &OutImageName, "Output file name [FORMAT:]path" );
    cl.AddOption( Key( 'I', "input-warp" ), &InWarpName, "Input file name for average warp" );
    cl.AddOption( Key( 'O', "output-warp" ), &OutWarpName, "Output file name for average warp" );

    cl.AddSwitch( Key( 'A', "with-affine" ), &WriteIncludeAffine, true, "Include affine component of first warp in output" );

    cl.AddSwitch( Key( "nn" ), &Interpolation, cmtk::Interpolators::NEAREST_NEIGHBOR, "Use nearest neigbor interpolation" );
    cl.AddSwitch( Key( "cubic" ), &Interpolation, cmtk::Interpolators::CUBIC, "Use tri-cubic interpolation" );

    cl.AddCallback( Key( "replace-from" ), AddReplaceFrom, "Replace from pattern" );
    cl.AddCallback( Key( "replace-to" ), AddReplaceTo, "Replace to pattern" );

    cl.Parse( argc, const_cast<const char**>( argv ) );

    inFileList.push_back( cl.GetNext() );
    
    const char* next = cl.GetNextOptional();
    while ( next ) 
      {
      inFileList.push_back( next );
      next = cl.GetNextOptional();
      }
    }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }
  
  std::vector<cmtk::SplineWarpXform::SmartPtr> warpList;
  std::list<cmtk::SplineWarpXform::SmartPtr> splineWarpList;
  std::vector<cmtk::UniformVolume::SmartPtr> volumeList;
  cmtk::TypedStreamStudylist studylist;

  char* referenceStudy = NULL;
  size_t idx = 0;

  for ( std::list<const char*>::const_iterator inFile = inFileList.begin(); inFile != inFileList.end(); ++inFile ) 
    {
    if ( Verbose ) 
      cmtk::StdErr.printf( "Opening studylist %s.\n", *inFile );

    if ( !studylist.Read( *inFile ) )
      {
      cmtk::StdErr << "Unable to read studylist " << *inFile << "\n";
      throw cmtk::ExitException( 1 );
      }
    
    if ( ! idx ) 
      {
      referenceStudy = strdup( studylist.GetReferenceStudyPath() );
      } 
    else
      {
      if ( strcmp( studylist.GetReferenceStudyPath(), referenceStudy ) ) 
	{
	cmtk::StdErr.printf( "ERROR: Studylist #%d has a different reference study.\n", idx );
	//continue;
	throw cmtk::ExitException( 1 );
	}
      }
    
    cmtk::SplineWarpXform::SmartPtr splineWarp = cmtk::SplineWarpXform::SmartPtr::DynamicCastFrom( studylist.GetWarpXform() );
    if ( splineWarp ) 
      {
      std::string actualPath = cmtk::StrReplace( studylist.GetFloatingStudyPath(), ReplaceMap );
      
      cmtk::UniformVolume::SmartPtr nextVolume;

      if ( OutImageName && !Jacobian ) 
	// no out image, no need for input images
	nextVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( actualPath.c_str(), Verbose ) );

      if ( OutImageName && !Jacobian && !nextVolume ) 
	{
	cmtk::StdErr.printf( "WARNING: Cannot read volume %s in studylist #%d.\n", actualPath.c_str(), idx );
	} 
      else
	{	
	// everything worked out, so let's push this item onto lists
	if ( SetPaddingValue )
	  {
	  nextVolume->GetData()->SetPaddingValue( PaddingValue );
	  if ( ReplacePaddingValue )
	    {
	    nextVolume->GetData()->ReplacePaddingData( NewPaddingValue );
	    }
	  }

	warpList.push_back( splineWarp );
	splineWarpList.push_back( splineWarp );
	if ( OutImageName ) 
	  volumeList.push_back( nextVolume );
	// and increase valid warp counter.
	++idx;
	}
      } 
    else 
      {
      cmtk::StdErr.printf( "ERROR: Studylist #%d has no spline warp.\n", idx );
      throw cmtk::ExitException( 1 );
      }
    }
  
  std::string actualPath = cmtk::StrReplace( studylist.GetReferenceStudyPath(), ReplaceMap );
  cmtk::UniformVolume::SmartPtr refVolume( cmtk::VolumeIO::ReadOriented( actualPath.c_str(), Verbose ) );
  
  if ( ! refVolume ) 
    {
    cmtk::StdErr.printf( "ERROR: Cannot read reference volume %s.\n", actualPath.c_str() );
    throw cmtk::ExitException( 1 );
    }
  
  // everything worked out, so let's push this item onto lists
  if ( SetPaddingValue )
    {
    refVolume->GetData()->SetPaddingValue( PaddingValue );
    if ( ReplacePaddingValue )
      {
      refVolume->GetData()->ReplacePaddingData( NewPaddingValue );
      }
    }
  
  if ( AutoScale )
    {
    if ( Jacobian ) 
      {
      cmtk::StdErr << "Warning: we don't really use image data for Jacobian map -- \nSkipping auto normalization step.\n";
      }
    else 
      {
      cmtk::Types::DataItem meanRef, varRef;
      refVolume->GetData()->GetStatistics( meanRef, varRef );
      
      for ( size_t volIdx = 0; volIdx < volumeList.size(); ++volIdx ) 
	{
	cmtk::TypedArray::SmartPtr data( volumeList[volIdx]->GetData() );
	cmtk::Types::DataItem meanFlt, varFlt;
	data->GetStatistics( meanFlt, varFlt );
	
	if ( Verbose ) 
	  {
	  cmtk::StdErr.printf( "Before auto-rescale: [1] %f +/- %f, [2] %f +/- %f\n", meanRef, sqrt( varRef ), meanFlt, sqrt( varFlt ) );
	  }

	// auto rescale, that is, determine scaling factor and offset so that 
	// after scaling, the intensities in both images have the same mean and
	// standard deviation. Note that due to the squaring of values in 
	// std.dev. computation, the spread will not be exactly identical.
	const cmtk::Types::DataItem factor = sqrt(varRef) / sqrt(varFlt);
	data->Rescale( factor, meanRef - factor * meanFlt );
	
	if ( Verbose ) 
	  {
	  data->GetStatistics( meanFlt, varFlt );
	  cmtk::StdErr.printf( "After auto-rescale: [1] %f +/- %f, [2] %f +/- %f\n", meanRef, sqrt( varRef ), meanFlt, sqrt( varFlt ) );
	  }
	}
      }
    }
  
  if ( Labels ) 
    {
    refVolume->GetData()->SetDataClass( cmtk::DATACLASS_LABEL );
    }

  cmtk::WarpXform::SmartPtr warpXform;
  if ( InWarpName ) 
    {
    cmtk::ClassStream stream( InWarpName, cmtk::ClassStream::READ );
    stream >> warpXform;
    } 
  else
    {
    cmtk::SmartPointer<cmtk::SplineActiveDeformationModel> adm( new cmtk::SplineActiveDeformationModel( splineWarpList, 0, IncludeScale, IncludeReferenceModel ) );
    adm->Compose();  
    warpXform = cmtk::WarpXform::SmartPtr::DynamicCastFrom( adm );
    }
  
#ifdef CMTK_BUILD_MPI
  if ( mpiRank == 0 )
#endif
    {
    if ( OutWarpName ) 
      {
      if ( WriteIncludeAffine ) 
	{
	warpXform->ReplaceInitialAffine( (*(warpList.begin()))->GetInitialAffineXform() );
	}
      
      cmtk::ClassStream stream( OutWarpName, cmtk::ClassStream::WRITE );
      stream << warpXform;
      }
    }

  if ( OutImageName ) 
    {
    cmtk::ProgressConsole progressIndicator;
    
    cmtk::ReformatVolume reformat;
    reformat.SetReferenceVolume( refVolume );
    reformat.SetWarpXform( warpXform );
    reformat.SetInterpolation( Interpolation );
    reformat.SetUsePaddingValue( HavePadOutValue );    
    reformat.SetPaddingValue( PadOutValue );
    if ( UserDataType != cmtk::TYPE_NONE )
      reformat.SetUserDataType( UserDataType );
    
    cmtk::UniformVolume::SmartPtr average;
    if ( Jacobian ) 
      {
      average = cmtk::UniformVolume::SmartPtr( reformat.GetTransformedReferenceJacobianAvg( &warpList, NULL /*origin*/, IncludeReferenceData ) );
      } 
    else 
      {
      average = cmtk::UniformVolume::SmartPtr( reformat.GetTransformedReference( &warpList, &volumeList, NULL /*origin*/, IncludeReferenceData ) );
      }
    
#ifdef CMTK_BUILD_MPI
    if ( mpiRank == 0 )
#endif
      cmtk::VolumeIO::Write( *average, OutImageName, Verbose );
    }
  
  free( referenceStudy );
  
#ifdef CMTK_BUILD_MPI
  MPI::Finalize();
#endif

  return 0;
}

int
main( int argc, char* argv[] )
{
  int exitCode = 0;
  try
    {
    exitCode = doMain( argc, argv );
    }
  catch ( const cmtk::ExitException& ex )
    {
    exitCode = ex.ExitCode();
    }
  return exitCode;
}
