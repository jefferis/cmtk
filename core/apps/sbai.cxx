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
#include <System/cmtkStrUtility.h>
#include <System/cmtkTimers.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkXformUniformVolume.h>
#include <Base/cmtkAffineXformUniformVolume.h>
#include <Base/cmtkSplineWarpXformUniformVolume.h>

#include <IO/cmtkFileFormat.h>
#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkTypedStreamStudylist.h>

#include <Segmentation/cmtkLabelCombinationShapeBasedAveragingInterpolation.h>

#include <math.h>
#include <vector>

#ifdef CMTK_USE_GCD
#  include <dispatch/dispatch.h>
#endif

const char *OutputFileName = "sbai.nii";

std::vector<std::string> InputFileVector;

unsigned short NumberOfLabels = 0;

const char* ReplaceFrom;
std::map<std::string,std::string> ReplaceMap;

bool PaddingFlag = false;
float PaddingValue = 0;

void AddReplaceFrom( const char* arg )
{
  ReplaceFrom = arg;
}

void AddReplaceTo( const char* arg )
{
  ReplaceMap[std::string(ReplaceFrom)] = std::string( arg );
}

void
AddVolumeStudyVector
( const char* listName, std::vector<cmtk::UniformVolume::SmartConstPtr>& volumeVector,
  std::vector<cmtk::XformUniformVolume::SmartConstPtr>& xformVector,
  cmtk::UniformVolume::SmartConstPtr& referenceVolume )
{
  cmtk::DebugOutput( 1 ) << "Opening studylist " << listName << ".\n";

  cmtk::TypedStreamStudylist studylist;
  studylist.Read( listName );

  if ( ! referenceVolume )
    {
    const std::string actualPath = cmtk::StrReplaceByRules( studylist.GetReferenceStudyPath(), ReplaceMap );
    referenceVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadGridOriented( actualPath ) );
    if ( ! referenceVolume )
      {
      cmtk::StdErr << "WARNING: could not read reference volume " << actualPath.c_str() << "\n";
      return;
      }
    }

  const std::string actualPath = cmtk::StrReplaceByRules( studylist.GetFloatingStudyPath(), ReplaceMap );
  cmtk::UniformVolume::SmartPtr floatingVolume( cmtk::VolumeIO::ReadOriented( actualPath ) );

  if ( PaddingFlag )
    {
    floatingVolume->GetData()->SetPaddingValue( PaddingValue );
    }

  volumeVector.push_back( floatingVolume );

  cmtk::SplineWarpXform::SmartConstPtr splinexform( cmtk::SplineWarpXform::SmartConstPtr::DynamicCastFrom( studylist.GetWarpXform() ) );
  if ( !splinexform )
    {
    cmtk::AffineXform::SmartConstPtr affinexform( studylist.GetAffineXform() );

    if ( !affinexform )
      {
      cmtk::StdErr << "WARNING: no transformation in studylist " << listName << "\n";
      return;
      }
    else
      {
      xformVector.push_back( cmtk::XformUniformVolume::SmartConstPtr( new cmtk::AffineXformUniformVolume( *(referenceVolume), *(affinexform) ) ) );
      }
    }
  else
    {
    xformVector.push_back( cmtk::XformUniformVolume::SmartConstPtr( new cmtk::SplineWarpXformUniformVolume( *(referenceVolume), splinexform ) ) );
    }
}

int
doMain ( const int argc, const char* argv[] ) 
{
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Shape-Based Averaging and interpolation of label images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Average segmentations (label fields) using the Euclidean Distance Transform. This tool performs joint interpolation and averaging by interpolating from the EDT. "
		       "This requires that the inputs are transformations from the same fixed to (not necessarily) different moving images. EDT computation is done in the space of each moving image. "
		       "See http://dx.doi.org/10.1109/TIP.2006.884936 for details of the underlying algorithm." );
    
    cl.AddParameterVector( &InputFileVector, "InputImageVector", "Input transformation file names." );
    
    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Input", "Input Options" );
    cl.AddOption( Key( 'n', "num-labels" ), &NumberOfLabels, "Number of labels. It is assumed that only values [0..num] occur in the images" );
    cl.AddOption( Key( 'p', "padding" ), &PaddingValue, "Padding value in input image", &PaddingFlag );
    
    cl.AddCallback( Key( "replace-from" ), AddReplaceFrom, "Replace from pattern" );
    cl.AddCallback( Key( "replace-to" ), AddReplaceTo, "Replace to pattern" );
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
  cmtk::UniformVolume::SmartConstPtr referenceVolume;

  std::vector<cmtk::XformUniformVolume::SmartConstPtr> xformVector;
  for ( std::vector<std::string>::const_iterator it = InputFileVector.begin(); it != InputFileVector.end(); ++it ) 
    {
    AddVolumeStudyVector( it->c_str(), volumeVector, xformVector, referenceVolume );
    }

  const double timeBaseline = cmtk::Timers::GetTimeProcess();
  cmtk::TypedArray::SmartPtr avgArray = cmtk::TypedArray::SmartPtr( cmtk::LabelCombinationShapeBasedAveragingInterpolation( volumeVector, xformVector, referenceVolume, NumberOfLabels ).GetResult() );
  cmtk::DebugOutput( 1 ).GetStream().printf( "Time %f sec\n", cmtk::Timers::GetTimeProcess() - timeBaseline );

  cmtk::UniformVolume::SmartPtr volume = referenceVolume->CloneGrid();  
  volume->SetData( avgArray );
  cmtk::VolumeIO::Write( *volume, OutputFileName );

  return 0;
}

#include "cmtkSafeMain"
