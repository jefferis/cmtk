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

#include <cmtkValueSequence.h>
#include <cmtkTemplateArray.h>
#include <cmtkSplineWarpXform.h>

#include <cmtkTypedStreamStudylist.h>
#include <cmtkVolumeIO.h>

#include <stdio.h>
#include <string.h>

#ifdef HAVE_LIMITS_H
#  include <limits.h>
#endif

#include <iostream>
#include <list>

bool Verbose = false;

bool WarpOnly = false;

int Direction = -1;
int ScanDirection = 0;

bool UseMask = false;
const char *MaskFile = NULL;
bool MultiValueMask = false;

bool UseThreshold = false;
float Threshold = 0;
const char* OutImagePath = NULL;

std::list<const char*> InListNames;

class AxesReorder 
{
public:
  virtual ~AxesReorder() {}
  virtual int NewX( const int x, const int y, const int z ) const = 0;
  virtual int NewY( const int x, const int y, const int z ) const = 0;
  virtual int NewZ( const int x, const int y, const int z ) const = 0;

protected:
  int DimsX;
  int DimsY;
  int DimsZ;
};

class AxesReorderAxial : public AxesReorder 
{
public:
  virtual ~AxesReorderAxial() {}
  AxesReorderAxial( const int dims[3] ) {
    DimsX = dims[0]-1; DimsY = dims[1]-1; DimsZ = dims[2]-1;
  }

  virtual int NewX( const int x, const int, const int ) const { return x; }
  virtual int NewY( const int, const int y, const int ) const { return y; }
  virtual int NewZ( const int, const int, const int z ) const { return z; }
};

class AxesReorderSagittal : public AxesReorder 
{
public:
  virtual ~AxesReorderSagittal() {}
  AxesReorderSagittal( const int dims[3] ) {
    DimsX = dims[0]-1; DimsY = dims[1]-1; DimsZ = dims[2]-1;
  }

  virtual int NewX( const int, const int, const int z ) const { return z; }
  virtual int NewY( const int x, const int, const int ) const { return x; }
  virtual int NewZ( const int, const int y, const int ) const { return DimsY-y; }
};

class AxesReorderCoronal : public AxesReorder 
{
public:
  virtual ~AxesReorderCoronal() {}
  AxesReorderCoronal( const int dims[3] ) {
    DimsX = dims[0]-1; DimsY = dims[1]-1; DimsZ = dims[2]-1;
  }

  virtual int NewX( const int x, const int, const int ) const { return x; }
  virtual int NewY( const int, const int, const int z ) const { return DimsZ-z; }
  virtual int NewZ( const int, const int y, const int ) const { return y; }
};

int
main ( const int argc, const char* argv[] ) 
{
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Statistics on deformation fields" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] input0.list [input1.list ...]" );

    typedef cmtk::CommandLine::Key Key;    
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode on." );

    cl.AddOption( Key( 'm', "mask" ), &MaskFile, "Use mask file." );
    cl.AddOption( Key( 'M', "multi-mask" ), &MaskFile, "Use multi-valued mask file.", &MultiValueMask );

    cl.AddOption( Key( 't', "thresh" ), &Threshold, "Use threshold value.", &UseThreshold );

    cl.AddSwitch( Key( 'x' ), &Direction, 0, "Output x component" );
    cl.AddSwitch( Key( 'y' ), &Direction, 1, "Output y component" );
    cl.AddSwitch( Key( 'z' ), &Direction, 2, "Output z component" );
    cl.AddOption( Key( 'o', "output-file" ), &OutImagePath, "Output image path [default: no image output]. "
		  "Use %d printf field for index when using multiple input transformations." );

    cl.AddSwitch( Key( 'j', "jacobian" ), &Direction, -2, "Jacobian determinant output." );
    cl.AddSwitch( Key( 'J', "jacobian-normalized" ), &Direction, -3, "Jacobian determinant output with global scaling correction." );

    cl.AddSwitch( Key( "axial" ), &ScanDirection, 0, "Axial input slice orientation." );
    cl.AddSwitch( Key( "sagittal" ), &ScanDirection, 1, "Sagittal input slice orientation." );
    cl.AddSwitch( Key( "coronal" ), &ScanDirection, 2, "Coronal input slice orientation." );

    cl.AddSwitch( Key( 'w', "warp-only" ), &WarpOnly, true, "Output warp component only (excluding affine)." );
    
    cl.Parse();

    const char* next = cl.GetNext();
    while ( next )
      {
      InListNames.push_back( next );
      next = cl.GetNextOptional();
      }
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  cmtk::UniformVolume::SmartPtr volume;

  cmtk::UniformVolume::SmartPtr maskVolume;
  bool maskValues[256];
  memset( maskValues, 0, sizeof( maskValues ) );

  if ( MaskFile && !maskVolume ) 
    {
    maskVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadOriented( MaskFile, Verbose ) );
    if ( ! maskVolume ) 
      {
      cmtk::StdErr << "Could not read mask volume " << MaskFile << "\n";
      exit(1);
      }
    
    if ( MultiValueMask ) 
      {
      if ( Verbose )
	cmtk::StdErr << "Setting up mask value table...\n";

      const cmtk::TypedArray* maskData = maskVolume->GetData();
      for ( size_t i = 0; i < maskData->GetDataSize(); ++i ) 
	{
	cmtk::Types::DataItem maskValue;
	if ( maskData->Get( maskValue, i ) ) 
	  {	  
	  maskValues[static_cast<byte>( maskValue )] = true;
	  }
	}
      }
    }

  int inListIdx = 0;
  cmtk::UniformVolume::SmartPtr fieldVolume;
  cmtk::ShortArray::SmartPtr fieldData;

  for ( std::list<const char*>::const_iterator it = InListNames.begin(); it != InListNames.end(); ++it ) 
    {
    if ( Verbose )
      fprintf( stderr, "Reading input studylist %s.\n", *it );

    cmtk::TypedStreamStudylist studylist( *it );
    const char* study0 = studylist.GetReferenceStudyPath();
    
    cmtk::WarpXform::SmartPtr warpXform = studylist.GetWarpXform();
    if ( ! warpXform ) continue;

    cmtk::AffineXform::SmartPtr initialAffine = warpXform->GetInitialAffineXform();
    
    if ( !volume )
      {
      volume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadGridOriented( study0, Verbose ) );
      if ( !volume ) 
	{
	cmtk::StdErr << "Could not read data volume " << study0 << "\n";
	exit( 1 );
	}

      fieldVolume = cmtk::UniformVolume::SmartPtr( volume->Clone( false /* copyData */ ) );
      fieldData = cmtk::ShortArray::SmartPtr( new cmtk::ShortArray( volume->GetNumberOfPixels() ) );
      fieldVolume->SetData( cmtk::TypedArray::SmartPtr::DynamicCastFrom( fieldData ) );
      }
    
    warpXform->RegisterVolume( volume );
    
    int planeAxis, xAxis, yAxis;
    AxesReorder *reorder;
    switch ( ScanDirection ) 
      {
      default:
      case 0:
	xAxis = cmtk::AXIS_X;
	yAxis = cmtk::AXIS_Y;
	planeAxis = cmtk::AXIS_Z;
	reorder = new AxesReorderAxial( volume->GetDims() );
	break;
      case 1:
	xAxis = cmtk::AXIS_Y;
	yAxis = cmtk::AXIS_Z;
	planeAxis = cmtk::AXIS_X;
	reorder = new AxesReorderSagittal( volume->GetDims() );
	break;
      case 2:
	xAxis = cmtk::AXIS_X;
	yAxis = cmtk::AXIS_Z;
	planeAxis = cmtk::AXIS_Y;
	reorder = new AxesReorderCoronal( volume->GetDims() );
	break;
      }
    
    fieldData->ClearArray();
    short *field = fieldData->GetDataPtrTemplate();
    
    for ( int maskIdx = 0; maskIdx < 256; ++maskIdx ) 
      {      
      if ( MultiValueMask || maskIdx ) 
	if ( ! maskValues[maskIdx] ) continue;
      
      if ( Verbose && MultiValueMask )
	cmtk::StdErr << "Processing mask value " << maskIdx << "...\n";
      
      cmtk::ValueSequence<double> normSeq, xSeq, ySeq, zSeq, JacobianSeq, EnergySeq;
      
      cmtk::Vector3D v0, v1;
      int offset = 0;
      for ( int plane=0; plane<volume->GetDims( planeAxis ); ++plane ) 
	{	
	if ( warpXform ) 
	  {
	  float globalScale = warpXform->GetGlobalScaling();
	  
	  for ( int y = 0; y < volume->GetDims( yAxis ); ++y )
	    for ( int x = 0; x < volume->GetDims( xAxis ); ++x, ++offset ) 
	      {
	      int newX = reorder->NewX(x,y,plane);
	      int newY = reorder->NewY(x,y,plane);
	      int newZ = reorder->NewZ(x,y,plane);
	      if ( UseThreshold ) 
		{
		cmtk::Types::DataItem data;
		if ( ! volume->GetDataAt( data, newX, newY, newZ ) ) 
		  continue;
		if ( data < Threshold ) continue;
		}
	      if ( UseMask ) 
		{
		cmtk::Types::DataItem data;
		if ( !maskVolume->GetDataAt( data, newX, newY, newZ ) ) 
		  continue;
		if ( MultiValueMask ) 
		  {
		  if ( static_cast<byte>( data ) != maskIdx ) continue;
		  } 
		else 
		  {
		  if ( data == 0 ) continue;
		  }
		}
	      
	      double jDet = 0;
	      warpXform->GetJacobianDeterminantSequence( &jDet, newX, newY, newZ );
	      if ( Direction == -3 ) jDet /= globalScale;
	      JacobianSeq.Proceed( jDet );

	      EnergySeq.Proceed( warpXform->GetGridEnergy( volume->GetGridLocation( newX, newY, newZ ) ) );
				 
	      if ( (Direction == -2)  || (Direction == -3) ) 
		{// Jacobian
		field[offset] = static_cast<short>( 100 * (jDet-1) );
		} 
	      else
		{
		volume->GetGridLocation( v0, newX, newY, newZ );
		if ( WarpOnly )
		  initialAffine->ApplyInPlaceNonVirtual( v0 );
		warpXform->GetTransformedGrid( v1, newX, newY, newZ );
		v1 -= v0;
		if ( Direction < 0 ) 
		  {
		  field[offset] = (short) ( 1000 * v1.EuclidNorm() );
		  } 
		else
		  {
		  field[offset] = (short) ( 1000 * fabs(v1.XYZ[Direction]) );
		  }
		normSeq.Proceed( v1.EuclidNorm() );
		xSeq.Proceed( v1.XYZ[0] );
		ySeq.Proceed( v1.XYZ[1] );
		zSeq.Proceed( v1.XYZ[2] );
		}
	      }
	  }
	}
      
      fprintf( stdout, "\n\nResults for studylist %s.\n", *it );
      fputs( "===============================================\n", stdout );
      
      if ( MultiValueMask )
	fprintf( stdout, "\nAnalysis for mask value %d:\n\n", maskIdx );
      
      fprintf( stdout, "GridEnergyIntegral %f\n"
	       "AvgGridEnergyIntegral %f\n"
	       "  variance %f, stddev %f, range %f %f\n", 
	       EnergySeq.GetSum(), EnergySeq.GetAverage(),
	       EnergySeq.GetVariance(), sqrt( EnergySeq.GetVariance() ),
	       EnergySeq.GetMinimum(), EnergySeq.GetMaximum() );

      fprintf( stdout, "JacobianDeterminantIntegral %f\n"
	       "AvgJacobianDeterminantIntegral %f\n"
	       "  variance %f, stddev %f, range %f %f\n", 
	       JacobianSeq.GetSum(), JacobianSeq.GetAverage(),
	       JacobianSeq.GetVariance(), sqrt( JacobianSeq.GetVariance() ),
	       JacobianSeq.GetMinimum(), JacobianSeq.GetMaximum() );

      fprintf( stdout, "\nInverseJacobianDeterminantIntegral %f\n"
	       "AvgInverseJacobianDeterminantIntegral %f\n"
	       "  variance %f, stddev %f, range %f %f\n", 
	       1.0 / JacobianSeq.GetSum(), 1.0 / JacobianSeq.GetAverage(),
	       1.0 / JacobianSeq.GetVariance(), 
	       1.0 / sqrt( JacobianSeq.GetVariance() ),
	       1.0 / JacobianSeq.GetMinimum(), 
	       1.0 / JacobianSeq.GetMaximum() );

      if ( Direction != -2 ) 
	{
	fprintf( stdout, "\nRegionVolumeVoxels %d\n"
		 "RegionVolumeCMM %f\n\n", xSeq.GetNValues(), 
		 xSeq.GetNValues() * volume->AverageVoxelVolume() );
	
	const char *fs = 
	  "%11s:  min = %6.3f mm  max = %6.3f mm "
	  "avg = %6.3f mm stddev = %6.3f mm\n"
	  "             |min| = %6.3f mm |max| = %6.3f mm, |delta|= %6.3f\n";
	fprintf( stdout, fs, "x-Direction", 
		 xSeq.GetMinimum(), xSeq.GetMaximum(), 
		 xSeq.GetAverage(), sqrt( xSeq.GetVariance()),
		 xSeq.GetMinimumAbs(), xSeq.GetMaximumAbs(),
		 xSeq.GetMaximum() - xSeq.GetMinimum() );
	fprintf( stdout, fs, "y-Direction", 
		 ySeq.GetMinimum(), ySeq.GetMaximum(), 
		 ySeq.GetAverage(), sqrt( ySeq.GetVariance()),
		 ySeq.GetMinimumAbs(), ySeq.GetMaximumAbs(),
		 ySeq.GetMaximum() - ySeq.GetMinimum() );
	fprintf( stdout, fs, "z-Direction", 
		 zSeq.GetMinimum(), zSeq.GetMaximum(), 
		 zSeq.GetAverage(), sqrt( zSeq.GetVariance()),
		 zSeq.GetMinimumAbs(), zSeq.GetMaximumAbs(),
		 zSeq.GetMaximum() - zSeq.GetMinimum() );
	fprintf( stdout, fs, "Total", 
		 normSeq.GetMinimum(), normSeq.GetMaximum(), 
		 normSeq.GetAverage(), sqrt( normSeq.GetVariance()),
		 normSeq.GetMinimumAbs(), normSeq.GetMaximumAbs(),
		 normSeq.GetMaximum() - normSeq.GetMinimum() );
	fputs( "\n", stdout );
	}
      }

    if ( OutImagePath ) 	 
      {
      char path[PATH_MAX];
      if ( snprintf( path, PATH_MAX, OutImagePath, inListIdx++ ) > PATH_MAX )
	{
	cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
 	}
      else
	{
	cmtk::VolumeIO::Write( fieldVolume, OutImagePath, Verbose );   
	}
      }
    }
}

