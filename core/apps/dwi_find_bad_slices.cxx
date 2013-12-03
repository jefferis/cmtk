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
#include <System/cmtkFileUtils.h>

#include <IO/cmtkVolumeIO.h>

#include <Base/cmtkRegionIndexIterator.h>

#include <vector>
#include <string>

int
doMain
( const int argc, const char *argv[] )
{
  std::vector<std::string> dwiImagePaths;

  int sliceAxis = cmtk::AXIS_Z;

  cmtk::Types::DataItem standardDeviations = 3.0;
  const char* outputDir = NULL;
  cmtk::Types::DataItem paddingValue = -1;
  cmtk::ScalarDataType outputDataType = cmtk::TYPE_NONE;

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

    cmtk::CommandLine::EnumGroup<cmtk::ScalarDataType>::SmartPtr typeGroup = 
      cl.AddEnum( "convert-to", &outputDataType, "Scalar data type for the output images. If your padding value is negative but your input data unsigned, for example, make sure to select a signed data type for the output. "
		  "By default, the output data type is the same as the input type.");
    typeGroup->AddSwitch( Key( "char" ), cmtk::TYPE_CHAR, "8 bits, signed" );
    typeGroup->AddSwitch( Key( "byte" ), cmtk::TYPE_BYTE, "8 bits, unsigned" );
    typeGroup->AddSwitch( Key( "short" ), cmtk::TYPE_SHORT, "16 bits, signed" );
    typeGroup->AddSwitch( Key( "ushort" ), cmtk::TYPE_USHORT, "16 bits, unsigned" );
    typeGroup->AddSwitch( Key( "int" ), cmtk::TYPE_INT, "32 bits signed" );
    typeGroup->AddSwitch( Key( "uint" ), cmtk::TYPE_UINT, "32 bits unsigned" );
    typeGroup->AddSwitch( Key( "float" ), cmtk::TYPE_FLOAT, "32 bits floating point" );
    typeGroup->AddSwitch( Key( "double" ), cmtk::TYPE_DOUBLE, "64 bits floating point\n" );

    cl.BeginGroup( "detection", "Bad Slice Detection" );
    cl.AddOption( Key( "stdev" ), &standardDeviations, "Threshold for bad slice identification in units of intensity standard deviations over all corresponding slices from the remaining diffusion images." );
    cl.EndGroup();

    cl.BeginGroup( "output", "Output Options" );
    cl.AddOption( Key( 'p', "padding-value" ), &paddingValue, "Padding value to replace data of detected bad slices in output images." );
    cl.AddOption( Key( 'o', "output-directory" ), &outputDir, "File system path for writing images with bad slices masked out (i.e., filled with a padding value)." );
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
  std::vector<cmtk::UniformVolume::SmartPtr> dwiImages( dwiImagePaths.size() );
  for ( size_t i = 0; i < dwiImagePaths.size(); ++i )
    {
    dwiImages[i] = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::Read( dwiImagePaths[i] ) );

    if ( i && ! dwiImages[0]->GridMatches( *dwiImages[i] ) )
      {
      cmtk::StdErr << "ERROR: geometry of image '" << dwiImagePaths[i] << "' does not match that of image '" << dwiImagePaths[0] << "'\n";
      throw cmtk::ExitException( 1 );
      }

    if ( outputDataType != cmtk::TYPE_NONE )
      {
      dwiImages[i]->SetData( dwiImages[i]->GetData()->Convert( outputDataType ) );
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

	// if we have an image output path given, mark bad slice using user-provided padding value
	if ( outputDir )
	  {
	  cmtk::DataGrid& image = *(dwiImages[i]);
	  const cmtk::DataGrid::RegionType sliceRegion = dwiImages[i]->GetSliceRegion( sliceAxis, slice );
	  for ( cmtk::RegionIndexIterator<cmtk::DataGrid::RegionType> it( sliceRegion ); it != it.end(); ++it )
	    {
	    image.SetDataAt( paddingValue, image.GetOffsetFromIndex( it.Index() ) );
	    }
	  }
	}
      }    
    }

  if ( outputDir )
    {
    for ( size_t i = 0; i < dwiImages.size(); ++i )
      {
      const std::string outputPath = std::string( outputDir ) + CMTK_PATH_SEPARATOR + cmtk::FileUtils::Basename( dwiImagePaths[i] );
      cmtk::VolumeIO::Write( *(dwiImages[i]), outputPath );
      }
    }

  return 0;
}

#include "cmtkSafeMain"
