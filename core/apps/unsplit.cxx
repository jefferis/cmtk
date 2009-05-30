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

#include <cmtkCommandLine.h>
#include <cmtkConsole.h>

#include <cmtkUniformVolume.h>
#include <cmtkVolumeIO.h>

#include <vector>
#include <list>

int
main( const int argc, const char* argv[] )
{
  bool verbose = false;

  std::list<const char*> inputFilePaths;
  const char* outputFilePath = NULL;
  int axis = 2;

  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Unsplit images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Join separate image stacks into a single interleaved image volume" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] inImage0 inImage1 ..." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &verbose, true, "Verbose operation" );

    cl.AddOption( Key( 'o', "output" ), &outputFilePath, "Path for output image." );

    cl.AddSwitch( Key( 'a', "axial" ), &axis, 2, "Interleaved axial images [default]" );
    cl.AddSwitch( Key( 'c', "coronal" ), &axis, 1, "Interleaved coronal images" );
    cl.AddSwitch( Key( 's', "sagittal" ), &axis, 0, "Interleaved sagittal images" );

    cl.AddSwitch( Key( 'x' ), &axis, 0, "Interleaved along x axis" );
    cl.AddSwitch( Key( 'y' ), &axis, 1, "Interleaved along y axis" );
    cl.AddSwitch( Key( 'z' ), &axis, 2, "Interleaved along z axis [default]" );

    cl.Parse();

    const char *next = cl.GetNext();
    while ( next )
      {
      inputFilePaths.push_back( next );
      next = cl.GetNextOptional();
      }
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  int stackDims[3] = { 0,0,0 };
  cmtk::Types::Coordinate stackDelta[3] = { 1,1,1 };
  std::vector<cmtk::UniformVolume::SmartPtr> volumes;
  for ( std::list<const char*>::const_iterator it = inputFilePaths.begin(); it != inputFilePaths.end(); ++it )
    {
    cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( *it, verbose ) );
    if ( ! volume || ! volume->GetData() )
      {
      cmtk::StdErr << "ERROR: Could not read image " << *it << "\n";
      return 1;
      }

    if ( volumes.size() )
      {
      // check image dimensions
      for ( int dim = 0; dim < 3; ++dim )
	{
	if ( dim == axis )
	  {
	  if ( (volume->Dims[dim] != volumes[0]->Dims[dim]) && (volume->Dims[dim]+1 != volumes[0]->Dims[dim]) )
	    {
	    cmtk::StdErr << "ERROR: interleaving dimension of image " << *it << " must be same as, or one smaller than first image's\n";
	    exit( 1 );
	    }
	  stackDims[dim] += volume->Dims[dim];
	  }
	else
	  {
	  if ( volume->Dims[dim] != volumes[0]->Dims[dim] )
	    {
	    cmtk::StdErr << "ERROR: in-plane dimensions of image " << *it << " do not match first image's\n";
	    exit( 1 );
	    }
	  }
	}
      }
    else
      // ! volumes.size() -> first image
      {
      // set dims and deltas; will modify later
      for ( int dim = 0; dim < 3; ++dim )
	{
	stackDims[dim] = volume->Dims[dim];
	stackDelta[dim] = volume->Delta[dim];
	}
      }
    volumes.push_back( volume );
    }

  stackDelta[axis] = volumes[0]->Delta[axis] / volumes.size();
  
  if ( verbose )
    {
    cmtk::StdErr << "Stacked image will have dimensions " << stackDims[0] << "x" << stackDims[1] << "x" << stackDims[2] << "\n";
    cmtk::StdErr << "Stacked image will have pixel size " << stackDelta[0] << "x" << stackDelta[1] << "x" << stackDelta[2] << "\n";
    }

  cmtk::UniformVolume::SmartPtr stacked( new cmtk::UniformVolume( stackDims, stackDelta[0], stackDelta[1], stackDelta[2] ) );
  stacked->CreateDataArray( volumes[0]->GetData()->GetType() );

  int toSlice = 0;
  for ( int fromSlice = 0; fromSlice < volumes[0]->Dims[axis]; ++fromSlice )
    {
    for ( size_t fromVolume = 0; fromVolume < volumes.size(); ++fromVolume, ++toSlice )
      {
      cmtk::ScalarImage::SmartPtr slice( volumes[fromVolume]->GetOrthoSlice( axis, fromSlice ) );
      stacked->SetOrthoSlice( axis, toSlice, slice );
      }
    }

  // get origin and orientation from first input image
  cmtk::AffineXform::MatrixType xformMatrix = volumes[0]->GetImageToPhysicalMatrix();
  // and copy to output
  stacked->m_IndexToPhysicalMatrix *= xformMatrix;
  stacked->m_MetaInformation = volumes[0]->m_MetaInformation;
  
  if ( outputFilePath )
    {
    cmtk::VolumeIO::Write( stacked, outputFilePath, verbose );
    }

  return 0;
}

