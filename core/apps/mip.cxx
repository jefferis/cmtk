/*
//
//  Copyright 2004-2010 SRI International
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <cmtkCommandLine.h>
#include <cmtkMathUtil.h>

#include <cmtkPGM.h>
#include <cmtkVolumeIO.h>
#include <cmtkDataGrid.h>
#include <cmtkAccumulatorMax.h>
#include <cmtkAccumulatorAvg.h>

#include <iostream>

bool Verbose = false;

double Black = cmtk::MathUtil::GetDoubleNaN();
double White = cmtk::MathUtil::GetDoubleNaN();

int Axis = 2;

bool Absolute = false;
bool AdjustAspect = false;

const char* InFileName = NULL;
const char* OutFileName = NULL;

bool SumProjection = false;
bool Write16Bit = false;

bool
ParseCommandLine( const int argc, const char* argv[] )
{
  try 
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Maximum Intensity Projection" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Compute 2D Maximum Intensity Projection image from a 3D volume image along one of the coordinate directions" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] input output" );
    
    typedef cmtk::CommandLine::Key Key;    
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode." );
    cl.AddSwitch( Key( 'a', "absolute" ), &Absolute, true, "Use absolute intensity values." );
    cl.AddSwitch( Key( 'A', "adjust" ), &AdjustAspect, true, "Adjust aspect ratio of output image (square pixels).\n" );
    
    cl.AddSwitch( Key( "axial" ), &Axis, 2, "\tAxial projection image." );
    cl.AddSwitch( Key( "coronal" ), &Axis, 1, "Coronal projection image." );
    cl.AddSwitch( Key( "sagittal" ), &Axis, 0, "Sagittal projection image.\n" );
    
    cl.AddOption( Key( 'b', "black" ), &Black, "Set black colormap value." );
    cl.AddOption( Key( 'w', "white" ), &White, "Set white colormap value." );
    
    cl.AddSwitch( Key( 'S', "sum" ), &SumProjection, true, "Sum projection." );
    cl.AddSwitch( Key( "write-16bit" ), &Write16Bit, true, "Write output as non-standard 16bit PGM image (doubles dynamic range but cannot be properly read by all readers)." );
    
    if ( ! cl.Parse() ) return false;
    
    InFileName = cl.GetNext();
    OutFileName = cl.GetNext();
    }
  catch ( cmtk::CommandLine::Exception& ex ) 
    {
    cmtk::StdErr << ex << "\n";
    return false;
    }
  
  return true;
}

int main ( const int argc, const char* argv[] ) 
{
  if ( ! ParseCommandLine( argc, argv ) ) 
    {
    return 1;
    }
  
  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( InFileName, Verbose ) );
  
  if ( Absolute ) 
    volume->GetData()->MakeAbsolute();
  
  cmtk::ScalarImage::SmartPtr mip;
  if ( SumProjection )
    mip = cmtk::ScalarImage::SmartPtr( volume->ComputeProjection< cmtk::Accumulators::Avg<cmtk::Types::DataItem> >( Axis ) );
  else
    mip = cmtk::ScalarImage::SmartPtr( volume->ComputeProjection< cmtk::Accumulators::Max<cmtk::Types::DataItem> >( Axis ) );
  if ( AdjustAspect ) 
    mip->AdjustAspect();
  
  cmtk::Types::DataItem min, max;
  mip->GetPixelData()->GetRange( min, max );
  if ( isnan( Black ) )
    Black = min;
  if ( isnan( White ) )
    White = max;

  if ( Write16Bit )
    cmtk::PGM::Write16bit( OutFileName, mip, Black, White );
  else
    cmtk::PGM::Write( OutFileName, mip, Black, White );
  
  return 0;
}

