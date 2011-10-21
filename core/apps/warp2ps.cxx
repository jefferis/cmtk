/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <iostream>

#include <Base/cmtkVector3D.h>
#include <Base/cmtkSplineWarpXform.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformIO.h>

const char* RefFileName = NULL;

const char* InListName = NULL;
bool InvertXform = false;

const char* MinusXformPath = NULL;
bool InvertMinusXform = false;

float ScaleFactor = 1.0;
int SamplingFactor = 1;
float LineWidth = 0;

int SliceIndex = -1;
int Axis = 2;

bool DrawBox = false;

bool Crop = false;
int CropRegion[4] = { 0, 0, 0, 0 };

void
SetCropRegion( const char* arg )
{
  if ( 4 == sscanf( arg, "%d,%d:%d,%d", &CropRegion[0], &CropRegion[1], &CropRegion[2], &CropRegion[3] ) )
    Crop = true;
  else
    Crop = false;
}

void DrawLine( const std::vector<double>& outputX, const std::vector<double>& outputY )
{
  if ( outputX.empty() )
    return;
  
  std::cout << outputX[0] << " " << outputY[0] << " m\n";
  for ( size_t i = 1; i < outputX.size(); ++i )
    std::cout << outputX[i] << " " << outputY[i] << " l\n";
  std::cout << "s\n";
}

int
doMain ( const int argc, const char* argv[] )
{
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Deformation to PostScript" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Write deformation field as deformed grid in PostScript format for visualization and illustration" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, warp2ps "[options] referenceImage xform" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( "invert-xform" ), &InvertXform, true, "Invert the main transformation" );
    cl.AddOption( Key( 'm', "minus-xform" ), &MinusXformPath, "Subtract this transformation (e.g., to exclude initial affine component)" );
    cl.AddSwitch( Key( "invert-minus-xform" ), &InvertMinusXform, true, "Invert the minus transformation" );

    cl.AddSwitch( Key( 'Z', "axis-z" ), &Axis, 2, "Slice orthogonal to z-axis (axial, default)" );
    cl.AddSwitch( Key( 'Y', "axis-y" ), &Axis, 1, "Slice orthogonal to y-axis (sagittal)" );
    cl.AddSwitch( Key( 'X', "axis-x" ), &Axis, 0, "Slice orthogonal to x-axis (coronal)" );
    cl.AddOption( Key( "slice" ), &SliceIndex, "Index of z-slice (default: center )");

    cl.AddOption( Key( "line-width" ), &LineWidth, "Line width [default: 0]" );
    cl.AddOption( Key( "sampling" ), &SamplingFactor, "Sampling factor (default: 1)");
    cl.AddOption( Key( "scale" ), &ScaleFactor, "Coordinage scale factor (default: 1.0)");
    cl.AddSwitch( Key( 'b', "box" ), &DrawBox, true, "Draw reference bounding box" );
    cl.AddCallback( Key( 'C', "crop" ), SetCropRegion, "Set in-plane crop region as 'x0,y0:x1,y1'" );

    if ( ! cl.Parse( argc, argv ) ) return 1;
    
    RefFileName = cl.GetNext();
    InListName = cl.GetNext();
    }
  catch ( const cmtk::CommandLine::Exception& ex ) 
    {
    cmtk::StdErr << ex;
    return 1;
    }
  

  cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( InListName ) );;
  if ( ! xform ) 
    {
    cmtk::StdErr << "ERROR: could not read transformation " << InListName << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( RefFileName ) );
  if ( ! volume ) 
    {
    cmtk::StdErr << "ERROR: could not read image " << RefFileName << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::Xform::SmartPtr minusXform( cmtk::Xform::SmartPtr::Null() );
  if ( MinusXformPath )
    {
    minusXform = cmtk::Xform::SmartPtr( cmtk::XformIO::Read( MinusXformPath ) );
    if ( ! minusXform ) 
      {
      cmtk::StdErr << "ERROR: could not read transformation " << MinusXformPath << "\n";
      throw cmtk::ExitException( 1 );
      }
    }

  int axisX = 0, axisY = 1, axisZ = 2;
  switch ( Axis )
    {
    default:
    case 0:
      axisX = 1;
      axisY = 2;
      axisZ = 0;
      break;
    case 1:
      axisX = 0;
      axisY = 2;
      axisZ = 1;
      break;
    case 2:
      axisX = 0;
      axisY = 1;
      axisZ = 2;
      break;
    }

  if ( ! Crop )
    {
    CropRegion[0] = 0;
    CropRegion[1] = 0;
    switch ( Axis )
      {
      default:
      case 0:
	CropRegion[2] = volume->GetDims()[cmtk::AXIS_Y];
	CropRegion[3] = volume->GetDims()[cmtk::AXIS_Z];
	break;
      case 1:
	CropRegion[2] = volume->GetDims()[cmtk::AXIS_X];
	CropRegion[3] = volume->GetDims()[cmtk::AXIS_Z];
	break;
      case 2:
	CropRegion[2] = volume->GetDims()[cmtk::AXIS_X];
	CropRegion[3] = volume->GetDims()[cmtk::AXIS_Y];
	break;
      }
    }

  cmtk::DebugOutput( 1 ).GetStream().printf( "Region: [%d,%d] ... [%d,%d]\n", CropRegion[0], CropRegion[1], CropRegion[2], CropRegion[3] );
  
  if ( SliceIndex < 0 )
    SliceIndex = volume->GetDims()[Axis] / 2;
  
  const cmtk::Types::Coordinate pixelSizeX = volume->Size[axisX] / (volume->GetDims()[axisX]-1);
  const cmtk::Types::Coordinate pixelSizeY = volume->Size[axisY] / (volume->GetDims()[axisY]-1);

  const cmtk::Types::Coordinate xmin = CropRegion[0] * pixelSizeX * ScaleFactor;
  const cmtk::Types::Coordinate ymin = CropRegion[1] * pixelSizeY * ScaleFactor;
  const cmtk::Types::Coordinate xmax = CropRegion[2] * pixelSizeX * ScaleFactor;
  const cmtk::Types::Coordinate ymax = CropRegion[3] * pixelSizeY * ScaleFactor;

  // PostScript header
  std::cout << "%!PS-Adobe-1.0\n";
  std::cout << "%%BoundingBox: " << xmin << " " << ymin << " " << xmax << " " << ymax << "\n";
  std::cout << "/m {moveto} def\n/l {lineto} def\n";
  std::cout << "/s {stroke} def\n/n {newpath} def\n";
  std::cout << "2 setlinejoin\n";
  std::cout << LineWidth << " setlinewidth\n";

  cmtk::UniformVolume::CoordinateVectorType v0, v1;
  std::vector<double> outputX;
  std::vector<double> outputY;

  const int SamplingFactorInLine = (SamplingFactor > 1) ? SamplingFactor / 2 : SamplingFactor;

  int idx[3] = { 0, 0, 0 };
  idx[axisZ] = SliceIndex;

  // vertical lines
  for ( idx[axisY] = CropRegion[1]; idx[axisY] < CropRegion[3]; idx[axisY] += SamplingFactor ) 
    {
    std::cout << "n\n";
    //      std::cout << "[0.7] 1 setdash\n";
    std::cout << LineWidth << " setlinewidth\n";
    
    outputX.resize( 0 );
    outputY.resize( 0 );
    for ( idx[axisX] = CropRegion[0]; idx[axisX] < CropRegion[2]; idx[axisX] += SamplingFactorInLine ) 
      {      
      bool success = true;
      v0 = v1 = volume->GetGridLocation( idx[0], idx[1], idx[2] );
      if ( xform )
	{
	if ( InvertXform )
	  success = success && xform->ApplyInverseInPlace( v1 );
	else
	  xform->ApplyInPlace( v1 );
	}
      if ( minusXform )
	{
	if ( InvertMinusXform )
	  success = success && minusXform->ApplyInverseInPlace( v0 );
	else
	  minusXform->ApplyInPlace( v0 );
	}

      if ( success )
	{
	v1 -= v0;
	v0 = volume->GetGridLocation( idx[0], idx[1], idx[2] );
	v1 += v0;
	
	outputX.push_back( ScaleFactor*v1[axisX] );
	outputY.push_back( ymax - ScaleFactor*v1[axisY] );
	}
      }
    DrawLine( outputX, outputY );
    }
  
  // horizontal lines
  for ( idx[axisX] = CropRegion[0]; idx[axisX] < CropRegion[2]; idx[axisX] += SamplingFactor ) 
    {
    std::cout << "n\n";
    //      std::cout << "[0.7] 1 setdash\n";
    std::cout << LineWidth << " setlinewidth\n";
    
    outputX.resize( 0 );
    outputY.resize( 0 );
    for ( idx[axisY] = CropRegion[1]; idx[axisY] < CropRegion[3]; idx[axisY] += SamplingFactorInLine ) 
      {
      bool success = true;
      v0 = v1 = volume->GetGridLocation( idx[0], idx[1], idx[2] );
      if ( xform )
	{
	if ( InvertXform )
	  success = success && xform->ApplyInverseInPlace( v1 );
	else
	  xform->ApplyInPlace( v1 );
	}
      if ( minusXform )
	{
	if ( InvertMinusXform )
	  success = success && minusXform->ApplyInverseInPlace( v0 );
	else
	  minusXform->ApplyInPlace( v0 );
	}

      if ( success )
	{
	v1 -= v0;
	v0 = volume->GetGridLocation( idx[0], idx[1], idx[2] );
	v1 += v0;
	
	outputX.push_back( ScaleFactor*v1[axisX] );
	outputY.push_back( ymax - ScaleFactor*v1[axisY] );
	}
      }
    DrawLine( outputX, outputY );
    }
  
  // reference box
  if ( DrawBox )
    {
    std::cout << LineWidth << " setlinewidth\n";
    std::cout << xmin << " " << ymin << " m "
	      << xmax << " " << ymin << " l "
	      << xmax << " " << ymax << " l "
	      << xmin << " " << ymax << " l "
	      << xmin << " " << ymin << " l\ns\n";
    }
  
// finalize
  std::cout << "showpage\n";

  return 0;
}

#include "cmtkSafeMain"
