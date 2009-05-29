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

#include <cmtkClassStream.h>
#include <cmtkClassStreamAffineXform.h>
#include <cmtkVolumeIO.h>

#include <cmtkSplineWarpXform.h>
#include <cmtkDeformationField.h>
#include <cmtkXformIO.h>
#include <cmtkXformList.h>

#include <stdio.h>

bool Verbose = false;
bool Mask = false;

const char* RefFileName = NULL;
const char *OutFileName = NULL;

const char* Downsample = NULL;

cmtk::Types::Coordinate InversionToleranceFactor = 0.1;

int
main ( const int argc, const char *argv[] ) 
{
  cmtk::UniformVolume::SmartPtr volume;

  cmtk::WarpXform::SmartPtr warpXform;
  cmtk::AffineXform::SmartPtr inverseInitial;

  cmtk::XformList xformList;

  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Transformation to Deformation Field" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Convert parametric rigid or nonrigid transformation to deformation field, sampled at pixel locations of a given refernece image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] outFile referenceImage [-i|--inverse] inList0 [[-i|--inverse] inList1 ...]" );
    
    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );
    cl.AddSwitch( Key( 'm', "mask" ), &Mask, true, "Use reference image pixels as a binary mask." );

    cl.AddOption( Key( "inversion-tolerance-factor" ), &InversionToleranceFactor, "Factor for numerical tolerance of B-spline inversion [multiples of minimum grid pixel size; default=0.1]" );
    cl.AddOption( Key( "downsample" ), &Downsample, "Downsample grid by factors 'x,y,z' or by single factor 'xyz'" );
 
    cl.Parse();

    OutFileName = cl.GetNext();
    RefFileName = cl.GetNext();
    const char* next = cl.GetNextOptional();
    while (next)
      {
      bool inverse = false;
      if ( !strcmp( next, "-i" ) || !strcmp( next, "--inverse" ) )
	{
	inverse = true;
	next = cl.GetNext();
	}

      cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( next, Verbose ) );
      xformList.Add( xform, inverse );
      next = cl.GetNextOptional();
      }
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e;
    exit( 1 );
    }

  if ( Mask )
    volume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( RefFileName, Verbose ) );
  else
    volume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadGridOriented( RefFileName, Verbose ) );
  if ( ! volume ) 
    {
    cmtk::StdErr << "Could not read reference volume " << RefFileName << "\n";
    exit(1);
    }          

  xformList.SetEpsilon( InversionToleranceFactor * volume->GetMinDelta() );
  
  if ( Downsample )
    {
    int factors[3] = { 1, 1, 1 };
    const size_t nFactors = sscanf( Downsample, "%d,%d,%d", factors, factors+1, factors+2 );
    if ( nFactors == 1 )
      {
      factors[1] = factors[2] = factors[0];
      }
    else
      {
      if ( nFactors != 3 )
	{
	cmtk::StdErr << "ERROR: downsampling factors must either be three integers, x,y,z, or a single integer\n";
	exit( 1 );
	}
      }
    volume = cmtk::UniformVolume::SmartPtr( volume->GetDownsampled( factors ) );
    }
  
  cmtk::DeformationField::SmartPtr dfield( new cmtk::DeformationField( volume ) );
  
  cmtk::Progress::SetTotalSteps( volume->GetDims( cmtk::AXIS_Z ) );
#pragma omp parallel for
  for ( int z = 0; z < volume->GetDims( cmtk::AXIS_Z ); ++z )
    {
    cmtk::Vector3D v0, v1;

    cmtk::Progress::SetProgress( z );

    size_t offset = 3 * z * volume->GetDims( cmtk::AXIS_X ) * volume->GetDims( cmtk::AXIS_Y );
    for ( int y = 0; y < volume->GetDims( cmtk::AXIS_Y ); ++y )
      {
      for ( int x = 0; x < volume->GetDims( cmtk::AXIS_X ); ++x, offset+=3 ) 
	{
	volume->GetGridLocation( v0, x, y, z );
	v1 = v0;

	bool invalid = true;
	if ( (!Mask) || (volume->GetDataAt( x, y, z ) > 0) )
	  {
	  invalid = !xformList.Apply( v1 );
	  }

	if ( !invalid )
	  v1 -= v0;
	else
	  v1.Set( 1e10, 1e10, 1e10 );
	
	dfield->m_Parameters[offset+0] = v1.XYZ[0];
	dfield->m_Parameters[offset+1] = v1.XYZ[1];
	dfield->m_Parameters[offset+2] = v1.XYZ[2];
	}
      }
    }
  cmtk::Progress::Done();

  cmtk::XformIO::Write( dfield, OutFileName, Verbose );
}

