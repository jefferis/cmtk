/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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
//  $Revision: 3844 $
//
//  $LastChangedDate: 2012-02-10 15:46:46 -0800 (Fri, 10 Feb 2012) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>

#include <IO/cmtkVolumeIO.h>

#include <vector>
#include <string>

int
doMain
( const int argc, const char *argv[] )
{
  std::vector<std::string> dwiImagePaths;

  int sliceAxis = cmtk::AXIS_Z;

  cmtk::Types::DataItem standardDeviations = 3.0;
  const char* outputPath = NULL;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Find bad slices in set of diffusion-weighted images." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool reads a set of 3D diffusion-weighted MR images and finds bad slices. A bad slice in a diffusion image is detected as one whose mean intensity is outside a specified "
		       "interval around the mean of the means of all corresponding slices from the remaining diffusion images." );
    
    cl.AddParameterVector( &dwiImagePaths, "DiffusionImagePaths", "List of file system paths for all diffusion-weighted images. This will usually not include any b=0 images, although it is possible to include these as well.");

    typedef cmtk::CommandLine::Key Key;

    cl.BeginGroup( "Input", "Input Options" );
    cmtk::CommandLine::EnumGroup<int>::SmartPtr sliceGroup = cl.AddEnum( "slice-orientation", &sliceAxis, "Define slice orientation of the diffusion images." );
    sliceGroup->AddSwitch( Key( "axial" ), (int)cmtk::AXIS_Z, "Axial slices" );
    sliceGroup->AddSwitch( Key( "sagittal" ),(int)cmtk::AXIS_X, "Sagittal slices" );
    sliceGroup->AddSwitch( Key( "coronal" ), (int)cmtk::AXIS_Y, "Coronal slices" );
    sliceGroup->AddSwitch( Key( "slice-x" ), (int)cmtk::AXIS_X, "X coordinate axis is slice direction" );
    sliceGroup->AddSwitch( Key( "slice-y" ), (int)cmtk::AXIS_Y, "Y coordinate axis is slice direction" );
    sliceGroup->AddSwitch( Key( "slice-z" ), (int)cmtk::AXIS_Z, "Z coordinate axis is slice direction" );
    cl.EndGroup();

    cl.BeginGroup( "detection", "Bad Slice Detection" );
    cl.AddOption( Key( "stdev" ), &standardDeviations, "Threshold for bad slice identification in units of intensity standard deviations over all corresponding slices from the remaining diffusion images." );
    cl.EndGroup();

    cl.BeginGroup( "output", "Output Options" );
    cl.AddOption( Key( 'o', "output" ), &outputPath, "File system path for the output file that contains the list of identified bad slices." );
    cl.EndGroup();

    cl.Parse( argc, argv );

    if ( dwiImagePaths.size() < 6 )
      {
      throw cmtk::CommandLine::Exception( "At least 6 diffusion images are needed for a DTI acquisition." );
      }
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }
  
  // read all diffusion images and make sure their grids match
  std::vector<cmtk::UniformVolume::SmartConstPtr> dwiImages( dwiImagePaths.size() );
  for ( size_t i = 0; i < dwiImagePaths.size(); ++i )
    {
    dwiImages[i] = cmtk::UniformVolume::SmartConstPtr( cmtk::VolumeIO::Read( dwiImagePaths[i] ) );

    if ( i && ! dwiImages[0]->GridMatches( *dwiImages[i] ) )
      {
      cmtk::StdErr << "ERROR: geometry of image '" << dwiImagePaths[i] << "' does not match that of image '" << dwiImagePaths[0] << "'\n";
      throw cmtk::ExitException( 1 );
      }
    }

  // loop over all slices
  for ( int slice = 0; slice < dwiImages[0]->m_Dims[sliceAxis]; ++slice )
    {
    // compute the mean intensity for this slice in each volume
    std::vector<cmtk::Types::DataItem> sliceMeans( dwiImages.size() );
    for ( size_t i = 0; i < dwiImages.size(); ++i )
      {
      cmtk::Types::DataItem variance = 0; // dummy variable; don't need variance
      dwiImages[i]->ExtractSlice( sliceAxis, slice )->GetData()->GetStatistics( sliceMeans[i], variance );
      }

    // test for each image whether this slice is bad in it
    for ( size_t i = 0; i < dwiImages.size(); ++i )
      {
      // get a vector of means without the current test image
      std::vector<cmtk::Types::DataItem> meansOtherImages;
      for ( size_t ii = 0; ii < dwiImages.size(); ++ii )
	{
	if ( i != ii )
	  {
	  meansOtherImages.push_back( sliceMeans[ii] );
	  }
	}

      const cmtk::Types::DataItem otherImagesMean = cmtk::MathUtil::Mean( meansOtherImages );
      const cmtk::Types::DataItem otherImagesSdev = sqrt( cmtk::MathUtil::Variance( meansOtherImages, otherImagesMean ) );

      const cmtk::Types::DataItem distance = fabs( sliceMeans[i] - otherImagesMean ) / otherImagesSdev;
      if ( distance > standardDeviations )
	{
	cmtk::DebugOutput( 2 ) << "Bad slice #" << slice << " in image #" << i << " mean=" << sliceMeans[i] << " distance=" << distance << " filename " << dwiImagePaths[i] << "\n";
	}
      }    
    }

  return 0;
}

#include "cmtkSafeMain"
