/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013 SRI International
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
#include <System/cmtkCommandLine.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkExitException.h>
#include <System/cmtkProgressConsole.h>

#include <IO/cmtkVolumeIO.h>

#include <Base/cmtkFilterVolume.h>
#include <Base/cmtkUnits.h>

#ifdef CMTK_BUILD_UNSTABLE
#  include <Unstable/cmtkFilterVolumeCoupe.h>
#endif

#ifdef CMTK_USE_SQLITE
#  include <Registration/cmtkImageXformDB.h>
#endif

#include <string.h>
#include <list>

bool Studholme = false;
bool Rohlfing = false;

#ifdef CMTK_BUILD_UNSTABLE
bool Coupe = false;
cmtk::Types::Coordinate CoupeBeta = 1;
//float CoupeBlockRadius = 0;  commented out while testing, this is defined in cmtk::FilterVolume.h
cmtk::Types::Coordinate CoupeWindowRadius = 5;
#endif // #ifdef CMTK_BUILD_UNSTABLE

bool Gaussian = false;
float GaussianWidth = 0.0;
float IntensityGaussianSigma = 0.0;

float Radius = 2.0;
float IntensityBinWidth = 10.0;

std::string InputFileName;
std::string OutputFileName;

std::string AverageFileName;
std::string SubjectFileName;
std::string MaskFileName;
std::list<std::string> ImageNameList;

#ifdef CMTK_USE_SQLITE
std::string updateDB;
#endif

int 
doMain( const int argc, const char* argv[] )
{
  cmtk::ProgressConsole progressIndicator( "Image Filtering" );
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Filter a volume image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool applies spatial filtering operators, including cnotent-sensitive opersators, based on selective Gaussian kernels." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "filter [options] input output\n"
		       "filter [options] [-s,--studholme] input output average subject img0 [img1...]\n"
		       "filter [options] [-R,--rohlfing] input output subject" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Filtering" );
    
    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( 'g', "gaussian" ), &GaussianWidth, "Gaussian filter width (sigma)", &Gaussian );
    cl.AddOption( Key( 'r', "radius" ), &Radius, "Filter radius (truncated outside)" );
    cl.AddOption( Key( 'm', "mask" ), &MaskFileName, "Binary mask file name" );

    cl.AddSwitch( Key( 's', "studholme" ), &Studholme, true, "Use Studholme's consistent filtering" );
    cl.AddOption( Key( 'b', "bin-width" ), &IntensityBinWidth, "Intensity bin width for consistent filtering" );

#ifdef CMTK_BUILD_UNSTABLE
    cl.AddSwitch( Key( 'C', "coupe" ), &Coupe, true, "Use Coupe's blockwise nonlocal means denoising filter" );
#endif
    cl.AddSwitch( Key( 'R', "rohlfing" ), &Rohlfing, true, "Use Rohlfing's single-image consistent filtering" );
    cl.AddOption( Key( 'G', "intensity-gaussian" ), &IntensityGaussianSigma, "Intensity gaussian sigma for consistent filtering" );

#ifdef CMTK_USE_SQLITE
    cl.BeginGroup( "Database", "Image/Transformation Database" );
    cl.AddOption( Key( "db" ), &updateDB, "Path to image/transformation database that should be updated with the newly created image." );
    cl.EndGroup();
#endif

    cl.Parse( argc, argv );

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
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }
  
  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( InputFileName ) );
  if ( !volume || !volume->GetData() )
    {
    cmtk::StdErr << "ERROR: Could not read volume " << InputFileName << "\n";
    throw cmtk::ExitException( 1 );
    }
  
  cmtk::TypedArray::SmartPtr maskData( NULL );
  if ( !MaskFileName.empty() ) 
    {
    cmtk::UniformVolume::SmartPtr maskVolume( cmtk::VolumeIO::ReadOriented( MaskFileName ) );
    if ( maskVolume )
      maskData = maskVolume->GetData();
    else
      {
      cmtk::StdErr << "ERROR: Could not read mask file " << MaskFileName << "\n";
      throw cmtk::ExitException( 1 );
      }
    }
  
  if ( Studholme ) 
    {
    cmtk::UniformVolume::SmartPtr average( cmtk::VolumeIO::ReadOriented( AverageFileName ) );
    if ( ! average || ! average->GetData() ) 
      {
      cmtk::StdErr	<< "ERROR: Could not read average anatomical file " << AverageFileName << "\n";
      throw cmtk::ExitException( 1 );
      }
    
    cmtk::UniformVolume::SmartPtr subject( cmtk::VolumeIO::ReadOriented( SubjectFileName ) );
    if ( ! subject || ! subject->GetData() ) 
      {
      cmtk::StdErr	<< "ERROR: Could not read subject anatomical file " << SubjectFileName << "\n";
      throw cmtk::ExitException( 1 );
      }

    std::list<cmtk::TypedArray::SmartPtr> imageList;
    for ( std::list<std::string>::const_iterator it = ImageNameList.begin(); it != ImageNameList.end(); ++it ) 
      {
      cmtk::UniformVolume::SmartPtr next( cmtk::VolumeIO::ReadOriented( *it ) );
      if ( ! next || ! next->GetData() ) 
	{
	cmtk::StdErr << "ERROR: Could not read subject anatomical file " << *it << "\n";
	throw cmtk::ExitException( 1 );
	}
      imageList.push_back( next->GetData() );
      }
    
    cmtk::TypedArray::SmartPtr filtered
      ( cmtk::FilterVolume::StudholmeFilter( volume,  subject->GetData(), average->GetData(), maskData, imageList, IntensityBinWidth, cmtk::Units::GaussianSigma( GaussianWidth ), Radius ) );
    volume->SetData( filtered );
    } 
  else
    {
    if ( Rohlfing ) 
      {
      cmtk::UniformVolume::SmartPtr subject( cmtk::VolumeIO::ReadOriented( SubjectFileName ) );
      if ( ! subject || ! subject->GetData() ) 
	{
	cmtk::StdErr << "ERROR: Could not read subject anatomical file " << SubjectFileName << "\n";
	throw cmtk::ExitException( 1 );
	}

      cmtk::TypedArray::SmartPtr filtered
	( cmtk::FilterVolume::RohlfingFilter( volume,  subject->GetData(), maskData, cmtk::Units::GaussianSigma( IntensityGaussianSigma ), cmtk::Units::GaussianSigma( GaussianWidth ), Radius ) );
      volume->SetData( filtered );
      }
    else
      {
      if ( Gaussian ) 
        {
        cmtk::TypedArray::SmartPtr filtered( cmtk::FilterVolume::GaussianFilter( volume, cmtk::Units::GaussianSigma( GaussianWidth ), Radius, maskData ) );
        volume->SetData( filtered );
        }
      else 
        {
#ifdef CMTK_BUILD_UNSTABLE
        if ( Coupe ) 
          {
          cmtk::TypedArray::SmartPtr filtered( cmtk::FilterVolumeCoupe::CoupeFilter( volume, static_cast<int>( CoupeWindowRadius ), CoupeBeta ) );
          volume->SetData( filtered );
          }
#endif // #ifdef CMTK_BUILD_UNSTABLE
        }
      }      
    }
  
  cmtk::VolumeIO::Write( *volume, OutputFileName );

#ifdef CMTK_USE_SQLITE
  if ( !updateDB.empty() )
    {
    cmtk::ImageXformDB db( updateDB );
    db.AddImage( OutputFileName, InputFileName  );
    }
#endif

  return 0;
}

#include "cmtkSafeMain"
