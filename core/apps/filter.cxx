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
#include <cmtkProgressConsole.h>

#include <string.h>
#include <list>

#include <cmtkVolumeIO.h>
#include <cmtkFilterVolume.h>

#ifdef CMTK_USE_SQLITE
#  include <cmtkImageXformDB.h>
#endif

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace filter
{
#endif
bool Verbose = false;

bool Studholme = false;
bool Rohlfing = false;

bool Coupe = false;
cmtk::Types::Coordinate CoupeBeta = 1;
//float CoupeBlockRadius = 0;  commented out while testing, this is defined in cmtk::FilterVolume.h
cmtk::Types::Coordinate CoupeWindowRadius = 5;

bool Gaussian = false;
float GaussianWidth = 0.0;
float IntensityGaussianSigma = 0.0;

float Radius = 2.0;
float IntensityBinWidth = 10.0;

const char* InputFileName = NULL;
const char* OutputFileName = NULL;

const char* AverageFileName = NULL;
const char* SubjectFileName = NULL;
const char* MaskFileName = NULL;
std::list<const char*> ImageNameList;

#ifdef CMTK_USE_SQLITE
const char* updateDB = NULL;
#endif

int main( int argc, char* argv[] )
{
  cmtk::ProgressConsole progressIndicator( "Image Filtering" );
  try 
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Filter a volume image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] input output\n"
		       "[options] [-s,--studholme] input output average subject img0 [img1...]\n"
		       "[options] [-R,--rohlfing] input output subject" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Filtering" );
    
    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddOption( Key( 'g', "gaussian" ), &GaussianWidth, "Gaussian filter width", &Gaussian );
    cl.AddOption( Key( 'r', "radius" ), &Radius, "Filter radius (truncated outside)" );
    cl.AddOption( Key( 'm', "mask" ), &MaskFileName, "Binary mask file name" );

    cl.AddSwitch( Key( 's', "studholme" ), &Studholme, true, "Use Studholme's consistent filtering" );
    cl.AddOption( Key( 'b', "bin-width" ), &IntensityBinWidth, "Intensity bin width for consistent filtering" );

    cl.AddSwitch( Key( 'C', "coupe" ), &Coupe, true, "Use Coupe's blockwise nonlocal means denoising filter" );
    cl.AddSwitch( Key( 'R', "rohlfing" ), &Rohlfing, true, "Use Rohlfing's single-image consistent filtering" );
    cl.AddOption( Key( 'G', "intensity-gaussian" ), &IntensityGaussianSigma, "Intensity gaussian sigma for consistent filtering" );

#ifdef CMTK_USE_SQLITE
    cl.BeginGroup( "Database", "Image/Transformation Database" );
    cl.AddOption( Key( "db" ), &updateDB, "Path to image/transformation database that should be updated with the newly created image." );
    cl.EndGroup();
#endif

    cl.Parse();

    InputFileName = cl.GetNext();
    OutputFileName = cl.GetNext();
    
    if ( Studholme ) 
      {
      AverageFileName = cl.GetNext();
      SubjectFileName = cl.GetNext();
      
      const char* next = cl.GetNext();
      while ( next ) 
	{
	ImageNameList.push_back( next );
	next = cl.GetNextOptional();
	}
      }
    else
      {
      if ( Rohlfing ) 
	{
	SubjectFileName = cl.GetNext();
	}
      }
    }
  catch ( cmtk::CommandLine::Exception e ) 
    {
    cmtk::StdErr << e;
    exit( 1 );
    }
  
  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( InputFileName, Verbose ) );
  if ( !volume || !volume->GetData() )
    {
    cmtk::StdErr << "ERROR: Could not read volume " << InputFileName << "\n";
    exit( 1 );
    }
  
  cmtk::TypedArray::SmartPtr maskData( NULL );
  if ( MaskFileName ) 
    {
    cmtk::UniformVolume::SmartPtr maskVolume( cmtk::VolumeIO::ReadOriented( MaskFileName, Verbose ) );
    if ( maskVolume )
      maskData = maskVolume->GetData();
    else
      {
      cmtk::StdErr << "ERROR: Could not read mask file " << MaskFileName << "\n";
      exit( 1 );
      }
    }
  
  if ( Studholme ) 
    {
    cmtk::UniformVolume::SmartPtr average( cmtk::VolumeIO::ReadOriented( AverageFileName, Verbose ) );
    if ( ! average || ! average->GetData() ) 
      {
      cmtk::StdErr	<< "ERROR: Could not read average anatomical file " << AverageFileName << "\n";
      exit( 1 );
      }
    
    cmtk::UniformVolume::SmartPtr subject( cmtk::VolumeIO::ReadOriented( SubjectFileName, Verbose ) );
    if ( ! subject || ! subject->GetData() ) 
      {
      cmtk::StdErr	<< "ERROR: Could not read subject anatomical file " << SubjectFileName << "\n";
      exit( 1 );
      }

    std::list<cmtk::TypedArray::SmartPtr> imageList;
    for ( std::list<const char*>::const_iterator it = ImageNameList.begin(); it != ImageNameList.end(); ++it ) 
      {
      cmtk::UniformVolume::SmartPtr next( cmtk::VolumeIO::ReadOriented( *it, Verbose ) );
      if ( ! next || ! next->GetData() ) 
	{
	cmtk::StdErr << "ERROR: Could not read subject anatomical file " << *it << "\n";
	exit( 1 );
	}
      imageList.push_back( next->GetData() );
      }
    
    cmtk::TypedArray::SmartPtr filtered
      ( cmtk::FilterVolume::StudholmeFilter( volume,  subject->GetData(), average->GetData(), maskData, imageList, IntensityBinWidth, GaussianWidth, Radius ) );
    volume->SetData( filtered );
    } 
  else
    {
    if ( Rohlfing ) 
      {
      cmtk::UniformVolume::SmartPtr subject( cmtk::VolumeIO::ReadOriented( SubjectFileName, Verbose ) );
      if ( ! subject || ! subject->GetData() ) 
	{
	cmtk::StdErr << "ERROR: Could not read subject anatomical file " << SubjectFileName << "\n";
	exit( 1 );
	}

      cmtk::TypedArray::SmartPtr filtered
	( cmtk::FilterVolume::RohlfingFilter( volume,  subject->GetData(), maskData, IntensityGaussianSigma, GaussianWidth, Radius ) );
      volume->SetData( filtered );
      }
    else
      {
      if ( Gaussian ) 
        {
        cmtk::TypedArray::SmartPtr filtered( cmtk::FilterVolume::GaussianFilter( volume, GaussianWidth, Radius, maskData ) );
        volume->SetData( filtered );
        }
      else 
        {
        if ( Coupe ) 
          {
          cmtk::TypedArray::SmartPtr filtered( cmtk::FilterVolume::CoupeFilter( volume, CoupeWindowRadius, CoupeBeta ) );
          volume->SetData( filtered );
          }
        }
      }      
    }
  
  cmtk::VolumeIO::Write( volume, OutputFileName, Verbose );

#ifdef CMTK_USE_SQLITE
  if ( updateDB )
    {
    cmtk::ImageXformDB db( updateDB );
    db.AddImage( OutputFileName, InputFileName  );
    }
#endif

}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace filter
} // namespace apps
} // namespace cmtk
#endif

