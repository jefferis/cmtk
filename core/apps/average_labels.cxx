/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkConsole.h>
#include <cmtkCommandLine.h>
#include <cmtkProgress.h>

#include <cmtkUniformVolume.h>
#include <cmtkVolumeIO.h>
#include <cmtkXformIO.h>

#include <queue>

bool Verbose = false;

typedef std::deque< std::pair<cmtk::Xform::SmartPtr,cmtk::UniformVolume::SmartPtr> > XVQueue;
XVQueue XformVolumeList;
cmtk::UniformVolume::SmartPtr ReferenceImage;

const char* OutputImagePath = "average_labels_output.hdr";

int
main( const int argc, const char* argv[] )
{
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Label image averaging" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Average co-registered label images using partial volumes" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] reference xform0 atlas0 [xform1 atlas1 ...]" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddOption( Key( 'o', "output" ), &OutputImagePath, "Output image path" );

    cl.Parse();
    
    const char* refImagePath = cl.GetNext();
    ReferenceImage = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadGridOriented( refImagePath, Verbose ) );
    if ( !ReferenceImage )
      {
      cmtk::StdErr << "ERROR: cannot read reference volume " << refImagePath << "\n";
      exit( 1 );
      }    

    const char* nextXform = cl.GetNext();
    const char* nextVolume = cl.GetNext();
    while ( nextXform && nextVolume )
      {
      cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( nextXform, Verbose ) );
      if ( !xform )
	{
	cmtk::StdErr << "ERROR: cannot read transformation " << nextXform << "\n";
	exit( 1 );
	}

      cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( nextVolume, Verbose ) );
      if ( !volume )
	{
	cmtk::StdErr << "ERROR: cannot read volume " << nextVolume << "\n";
	exit( 1 );
	}

      XformVolumeList.push_back( std::pair<cmtk::Xform::SmartPtr,cmtk::UniformVolume::SmartPtr>( xform, volume ) );

      nextXform = cl.GetNextOptional();
      nextVolume = cl.GetNextOptional();
      }
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  cmtk::TypedArray::SmartPtr result( cmtk::TypedArray::Create( cmtk::TYPE_SHORT, ReferenceImage->GetNumberOfPixels() ) );

  const int* dims = ReferenceImage->GetDims();
  cmtk::Progress::Begin( 0, dims[2], 1, "Label image averaging" );
#pragma omp parallel for
  for ( int z = 0; z < dims[2]; ++z )
    {
    cmtk::Progress::SetProgress( z );
    size_t offset = z * dims[0] * dims[1];
    float labelWeights[256];
    for ( int y = 0; y < dims[1]; ++y )
      for ( int x = 0; x < dims[0]; ++x, ++offset )
	{
	memset( labelWeights, 0, sizeof( labelWeights ) );
	cmtk::Vector3D v = ReferenceImage->GetGridLocation( x, y, z );

	for ( XVQueue::iterator it = XformVolumeList.begin(); it != XformVolumeList.end(); ++it )
	  {
	  const cmtk::Xform* xform = it->first;
	  const cmtk::UniformVolume* volume = it->second;

	  cmtk::Vector3D vx( v );
	  xform->ApplyInPlace( vx );
	  	  
	  cmtk::ProbeInfo probeInfo;
	  if ( volume->ProbeNoXform( probeInfo, vx ) )
	    for ( int corner = 0; corner < 8; ++corner )
	      labelWeights[static_cast<byte>( probeInfo.Values[corner] )] += probeInfo.GetWeight( corner );
	  }

	short maxLabel = -1;
	float maxLabelWeight = -1;

	for ( int l = 0; l < 256; ++l )
	  {
	  if ( labelWeights[l] > maxLabelWeight )
	    {
	    maxLabelWeight = labelWeights[l];
	    maxLabel = l;
	    }
	  else
	    {
	    if ( labelWeights[l] == maxLabelWeight )
	      {
	      maxLabel = -1;
	      }
	    }
	  }

	result->Set( maxLabel, offset );
	}  
    }
  cmtk::Progress::Done();

  ReferenceImage->SetData( result );

  cmtk::VolumeIO::Write( ReferenceImage, OutputImagePath, Verbose );
  return 0;
}

