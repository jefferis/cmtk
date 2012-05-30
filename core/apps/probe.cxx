/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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
//  $Revision: 3003 $
//
//  $LastChangedDate: 2011-03-17 11:18:47 -0700 (Thu, 17 Mar 2011) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkVector3D.h>
#include <Base/cmtkUniformVolumeInterpolatorBase.h>
#include <Base/cmtkUniformVolumeInterpolator.h>
#include <Base/cmtkUniformVolumeInterpolatorPartialVolume.h>
#include <Base/cmtkSincInterpolator.h>
#include <Base/cmtkLinearInterpolator.h>
#include <Base/cmtkCubicInterpolator.h>
#include <Base/cmtkNearestNeighborInterpolator.h>

#include <IO/cmtkVolumeIO.h>

#include <stdio.h>

cmtk::UniformVolumeInterpolatorBase*
MakeInterpolator( const cmtk::UniformVolume& volume, const cmtk::Interpolators::InterpolationEnum interpolation, const int interpolationWindowRadius )
{
  switch ( interpolation ) 
    {
    default:
    case cmtk::Interpolators::LINEAR:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear> TInterpolator;
    return new TInterpolator( volume );
    }
    case cmtk::Interpolators::CUBIC:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Cubic> TInterpolator;
    return new TInterpolator( volume );
    }
    case cmtk::Interpolators::NEAREST_NEIGHBOR:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::NearestNeighbor> TInterpolator;
    return new TInterpolator( volume );
    }
    case cmtk::Interpolators::COSINE_SINC:
    {
    switch ( interpolationWindowRadius )
      {
      case 2:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<2> > TInterpolator;
      return new TInterpolator( volume );  
      }
      case 3:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<3> > TInterpolator;
      return new TInterpolator( volume );  
      }
      case 4:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<4> > TInterpolator;
      return new TInterpolator( volume );  
      }
      case 5:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<5> > TInterpolator;
      return new TInterpolator( volume );  
      }      
      default:
	cmtk::StdErr.printf( "ERROR: Sinc window radius %d is not supported.\n", (int)interpolationWindowRadius );
      }
    }

    case cmtk::Interpolators::HAMMING_SINC:
    {
    switch ( interpolationWindowRadius )
      {
      case 2:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::HammingSinc<2> > TInterpolator;
      return new TInterpolator( volume );  
      }
      case 3:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::HammingSinc<3> > TInterpolator;
      return new TInterpolator( volume );  
      }
      case 4:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::HammingSinc<4> > TInterpolator;
      return new TInterpolator( volume );  
      }
      case 5:
      {
      typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::HammingSinc<5> > TInterpolator;
      return new TInterpolator( volume );  
      }
      default:
	cmtk::StdErr.printf( "ERROR: Sinc window radius %d is not supported.\n", (int)interpolationWindowRadius );
      }
    }

    case cmtk::Interpolators::PARTIALVOLUME:
    {
    typedef cmtk::UniformVolumeInterpolatorPartialVolume TInterpolator;
    return new TInterpolator( volume );
    }

    return NULL;
    }
}

typedef enum
{
  COORDINATES_INDEXED,
  COORDINATES_ABSOLUTE,
  COORDINATES_RELATIVE,
  COORDINATES_PHYSICAL
} CoordinateModeEnum;

int
doMain( const int argc, const char *argv[] )
{
  const char* inputImagePath = NULL;
  const char* readOrientation = "RAS";

  CoordinateModeEnum mode;

  cmtk::Interpolators::InterpolationEnum interpolation = cmtk::Interpolators::NEAREST_NEIGHBOR;
  int interpolationWindowRadius = 3;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Probe image data." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool prints pixel values or symbolic labels at a list of user-provided image coordinates." );

    typedef cmtk::CommandLine::Key Key;
    cmtk::CommandLine::EnumGroup<CoordinateModeEnum>::SmartPtr modeGroup = cl.AddEnum( "coordinates", &mode, "Coordinate specification mode." );
    modeGroup->AddSwitch( Key( "indexed" ), COORDINATES_INDEXED, "Use grid indexes to specify coordinates. For each dimension, the valid value range is [0,Dims-1]." );
    modeGroup->AddSwitch( Key( "absolute" ), COORDINATES_ABSOLUTE, "Use absolute volume coordinates. For each dimension, the valid range is [0,FOV]." );
    modeGroup->AddSwitch( Key( "relative" ), COORDINATES_RELATIVE, "Use relative volume coordinates. For each dimension, the valid range is [0,1]." );
    modeGroup->AddSwitch( Key( "physical" ), COORDINATES_PHYSICAL, "Use physical volume coordinates. "
			  "Each given location is transformed into image coordinates via the inverse of the images's index-to-physical space matrix." );

    cmtk::CommandLine::EnumGroup<cmtk::Interpolators::InterpolationEnum>::SmartPtr interpolationGroup = cl.AddEnum( "interpolation", &interpolation, "Image interpolation method." );
    interpolationGroup->AddSwitch( Key( "nn" ), cmtk::Interpolators::NEAREST_NEIGHBOR, "Nearest neighbor interpolation" );
    interpolationGroup->AddSwitch( Key( "linear" ), cmtk::Interpolators::LINEAR, "Trilinear interpolation" );
    interpolationGroup->AddSwitch( Key( "cubic" ), cmtk::Interpolators::CUBIC, "Tricubic interpolation" );
    interpolationGroup->AddSwitch( Key( "pv" ), cmtk::Interpolators::PARTIALVOLUME, "Partial volume interpolation" );
    interpolationGroup->AddSwitch( Key( "sinc-cosine" ), cmtk::Interpolators::COSINE_SINC, "Sinc interpolation with cosine window" );
    interpolationGroup->AddSwitch( Key( "sinc-hamming" ), cmtk::Interpolators::HAMMING_SINC, "Sinc interpolation with Hamming window" );

    cl.AddOption( Key( "sinc-window-radius" ), &interpolationWindowRadius, "Window radius for Sinc interpolation" );

    cl.AddSwitch( Key( "no-reorient" ), &readOrientation, static_cast<const char*>( NULL ), "Disable image reorientation into RAS alignment." );

    cl.AddParameter( &inputImagePath, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartPtr volume;
  if ( readOrientation )
    volume = cmtk::VolumeIO::ReadOriented( inputImagePath, readOrientation );
  else
    volume = cmtk::VolumeIO::Read( inputImagePath );
      
  if ( ! volume || ! volume->GetData() ) 
    {
    cmtk::StdErr << "ERROR: could not read image " << inputImagePath << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolumeInterpolatorBase::SmartPtr interpolator( MakeInterpolator( *volume, interpolation, interpolationWindowRadius ) );

  FILE *infile = stdin;
  FILE *outfile = stdout;

  while ( ! feof( infile ) )
    {
    double xyz[3];
    if ( 3 == fscanf( infile, "%lf %lf %lf", xyz, xyz+1, xyz+2 ) )
      {
      cmtk::UniformVolume::SpaceVectorType v = cmtk::UniformVolume::SpaceVectorType( xyz );
      
      switch ( mode )
	{
	case COORDINATES_INDEXED:
	  // absolute image coordinate is index times pixel size
	  v *= volume->m_Delta;
	  break;
	case COORDINATES_ABSOLUTE:
	  // nothing to do - lookup will be done by absolute coordinate
	  break;
	case COORDINATES_RELATIVE:
	  // absolute image coordinate is relative times volume size
	  v *= volume->m_Size;
	  break;
	case COORDINATES_PHYSICAL:
	  // absolute image coordinate is physical transformed by inverse image-to-physical matrix
	  v *= volume->GetImageToPhysicalMatrix().GetInverse();	  
	  break;
	}
      
      cmtk::Types::DataItem value;
      if ( interpolator->GetDataAt( v, value ) )
	{
	fprintf( outfile, "%lf\n", static_cast<double>( value ) );
	}
      else
	{
	fprintf( outfile, "NAN\n" );
	}
      }
    }

  return 0;
}

#include "cmtkSafeMain" 

