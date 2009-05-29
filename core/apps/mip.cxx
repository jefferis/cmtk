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

#include <cmtkCommandLine.h>
#include <cmtkMathUtil.h>

#include <cmtkPGM.h>
#include <cmtkVolumeIO.h>
#include <cmtkDataGrid.h>
#include <cmtkAccumulatorMax.h>
#include <cmtkAccumulatorAvg.h>

#include <iostream>

bool Verbose = false;

int Black = -1;
int White = -1;

int Axis = 2;

bool Absolute = false;
bool AdjustAspect = false;

const char* InFileName = NULL;
const char* OutFileName = NULL;

bool SumProjection = false;

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
  
  cmtk::Volume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( InFileName, Verbose ) );
  
  if ( Absolute ) volume->GetData()->MakeAbsolute();

  cmtk::ScalarImage::SmartPtr mip;
  if ( SumProjection )
    mip = cmtk::ScalarImage::SmartPtr( volume->ComputeProjection< cmtk::Accumulators::Avg<cmtk::Types::DataItem> >( Axis ) );
  else
    mip = cmtk::ScalarImage::SmartPtr( volume->ComputeProjection< cmtk::Accumulators::Max<cmtk::Types::DataItem> >( Axis ) );
  if ( AdjustAspect ) mip->AdjustAspect();
  
  cmtk::PGM::Write( OutFileName, mip, Black, White );
  
  return 0;
}

