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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkconfig.h>

#include <cmtkConsole.h>
#include <cmtkCommandLine.h>
#include <cmtkProgress.h>

#include <cmtkTypedArray.h>
#include <cmtkUniformVolume.h>

#include <cmtkOverlapMeasures.h>

#include <cmtkVolumeIO.h>

#include <vector>

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace overlap
{
#endif
bool Verbose = false;

bool PaddingInFlag = false;
cmtk::Types::DataItem PaddingInValue = 0;

std::vector<cmtk::TypedArray::SmartPtr> vectorOfDataArrays;
const char* MaskFileName = NULL;

unsigned int NumberOfLabels = 0;
unsigned int FirstLabel = 0;
bool ByLabel = false;

int
main ( int argc, char* argv[] ) 
{
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Overlap computation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Compute overlap measures between two or more images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] image0 image1 [image2 ...]" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddOption( Key( 'p', "pad-in" ), &PaddingInValue, "Define padding value for input images", &PaddingInFlag );
    cl.AddOption( Key( 'm', "mask" ), &MaskFileName, "Mask file for optional region-based analysis" );
    
    cl.AddOption( Key( 'N', "num-labels" ), &NumberOfLabels, "Number of label values [default: autodetect]" );
    cl.AddOption( Key( 'f', "first-label" ), &FirstLabel, "Index of first label value [default: 0]; labels below this one are ignored" );
    cl.AddSwitch( Key( "by-label" ), &ByLabel, true, "Analysis by label, i.e., separately for each label in given (detected) range" );
    
    if (!cl.Parse()) 
      {
      exit( 1 );
      }
        
    const char* next = cl.GetNext();
    while ( next )
      {
      cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( next, Verbose ) );
      if ( ! volume || ! volume->GetData() )
	{
	cmtk::StdErr << "ERROR: Could not read image " << next << "\n";
	exit( 1 );
	}

      if ( PaddingInFlag )
	{
	volume->GetData()->SetPaddingValue( PaddingInValue );
	}
      
      if ( vectorOfDataArrays.size() )
	{
	if ( vectorOfDataArrays[0]->GetDataSize() != volume->GetNumberOfPixels() )
	  {
	  cmtk::StdErr << "ERROR: pixel count of image " << next << " does not match other images.\n";
	  exit( 1 );
	  }
	}
      vectorOfDataArrays.push_back( volume->GetData() );
      
      next = cl.GetNextOptional();
      }
    }
  catch ( cmtk::CommandLine::Exception ex ) 
    {
    cmtk::StdErr << ex;
    return false;
    }

  cmtk::TypedArray::SmartPtr maskData;
  if ( MaskFileName ) 
    {
    cmtk::UniformVolume::SmartPtr maskVolume( cmtk::VolumeIO::ReadOriented( MaskFileName, Verbose ) );
    if ( ! maskVolume || ! maskVolume->GetData() )
      {
      cmtk::StdErr << "ERROR: Could not read mask image " << MaskFileName << "\n";
      exit( 1 );
      }
    maskData = maskVolume->GetData();
    }

  cmtk::ConsoleProgress progressIndicator;

  double overlapEqual = 0, overlapVolume = 0, overlapInverse = 0;
  
  const cmtk::OverlapMeasures overlapMeasures( vectorOfDataArrays );

  if ( !NumberOfLabels )
    {
    NumberOfLabels = 1 + overlapMeasures.GetMaxLabelValue() - FirstLabel;
    }

  if ( ByLabel )
    {
    fprintf( stdout, "Label\tO[eq]\tO[vol]\tO[invvol]\n" );
    for ( unsigned int label = FirstLabel; label < FirstLabel+NumberOfLabels; ++label )
      { 
      if ( overlapMeasures.ComputeGroupwiseOverlap( label, 1 /*numberOfLabels*/, overlapEqual, overlapVolume, overlapInverse ) > 0 )
	{
	fprintf( stdout, "%d\t%lf\t%lf\t %lf\n", label, overlapEqual, overlapVolume, overlapInverse  );
	}
      }
    }
  else
    {
    overlapMeasures.ComputeGroupwiseOverlap( FirstLabel, NumberOfLabels, overlapEqual, overlapVolume, overlapInverse );
    fprintf( stdout, "O[equal] = %lf\nO[volume] = %lf\nO[inverse_volume] = %lf\n", overlapEqual, overlapVolume, overlapInverse  );
    }

  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace overlap
} // namespace apps
} // namespace cmtk
#endif
