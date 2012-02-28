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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkconfig.h>

#include <System/cmtkConsole.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkTimers.h>

#include <Base/cmtkUniformVolume.h>

#include <IO/cmtkVolumeIO.h>

#include <Segmentation/cmtkLabelCombinationShapeBasedAveraging.h>

#include <vector>

#ifdef CMTK_USE_GCD
#  include <dispatch/dispatch.h>
#endif

const char *OutputFileName = "sba.nii";

std::vector<std::string> InputFileVector;

unsigned short NumberOfLabels = 0;

bool PaddingFlag = false;
float PaddingValue = 0;

void
AddVolumeFile
( const char* fileName, std::vector<cmtk::UniformVolume::SmartConstPtr>& volumeVector )
{
  cmtk::DebugOutput( 1 ) << "Opening image " << fileName << ".\n";
  
  cmtk::UniformVolume::SmartPtr nextVolume( cmtk::VolumeIO::ReadOriented( fileName ) );
  
  if ( PaddingFlag )
    {
    nextVolume->GetData()->SetPaddingValue( PaddingValue );
    }
  
  if ( nextVolume->GetData()->GetType() != cmtk::TYPE_USHORT )
    {
    cmtk::StdErr << "WARNING: converting data to 'unsigned short'\n";
    
    nextVolume->SetData( cmtk::TypedArray::SmartPtr( nextVolume->GetData()->Convert( cmtk::TYPE_USHORT ) ) );
    }
  volumeVector.push_back( nextVolume );
}

int
doMain ( const int argc, const char* argv[] ) 
{
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Shape-based Averaging of label images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Average segmentations (label fields) using the Euclidean Distance Transform. All input images must be in the same space. EDT is computed in this space also. "
		       "See http://dx.doi.org/10.1109/TIP.2006.884936 for details of the underlying algorithm." );

    cl.AddParameterVector( &InputFileVector, "InputImageVector", "Input image file names." );

    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Input", "Input Options" );
    cl.AddOption( Key( 'n', "num-labels" ), &NumberOfLabels, "Number of labels. It is assumed that only values [0..num] occur in the images" );
    cl.AddOption( Key( 'p', "padding" ), &PaddingValue, "Padding value in input image", &PaddingFlag );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( 'o', "output" ), &OutputFileName, "File name for output segmentation file." );
    cl.EndGroup();

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }

  std::vector<cmtk::UniformVolume::SmartConstPtr> volumeVector;
  for ( std::vector<std::string>::const_iterator it = InputFileVector.begin(); it != InputFileVector.end(); ++it ) 
    {
    AddVolumeFile( it->c_str(), volumeVector );
    }

  const double timeBaseline = cmtk::Timers::GetTimeProcess();

  cmtk::TypedArray::SmartPtr avgArray = cmtk::TypedArray::SmartPtr( cmtk::LabelCombinationShapeBasedAveraging( volumeVector, NumberOfLabels ).GetResult() );
  cmtk::DebugOutput( 1 ).GetStream().printf( "Time %f sec\n", cmtk::Timers::GetTimeProcess() - timeBaseline );
  
  cmtk::UniformVolume::SmartPtr volume = volumeVector[0]->CloneGrid();
  volume->SetData( avgArray );
  cmtk::VolumeIO::Write( *volume, OutputFileName );

  return 0;
}

#include "cmtkSafeMain"
