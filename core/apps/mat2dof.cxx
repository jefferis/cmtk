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
#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkVector.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkMatrix4x4.h>

#include <IO/cmtkStudyList.h>
#include <IO/cmtkClassStreamOutput.h>
#include <IO/cmtkClassStreamAffineXform.h>
#include <IO/cmtkClassStreamStudyList.h>

#include <fstream>
#include <iostream>

void
ReadMatrix( cmtk::Types::Coordinate (&matrix)[4][4], std::istream& stream, const bool matrix3x3 )
{
  memset( matrix, 0 , sizeof( matrix ) );

  if ( matrix3x3 )
    {
    for ( int row = 0; row < 3; ++row )
      for ( int col = 0; col < 3; ++col )
	stream >> matrix[col][row];
    matrix[3][3] = 1.0;
    }
  else
    {
    for ( int row = 0; row < 4; ++row )
      for ( int col = 0; col < 4; ++col )
	stream >> matrix[col][row];
    }
}

int
doMain( const int argc, const char* argv[] )
{
  const char* CenterStr = NULL;
  const char* OffsetStr = NULL;
  const char* TranslateStr = NULL;
  const char* PixelSizeStr = NULL;
  
  const char* Reference = "reference";
  const char* Floating = "floating";
  
  const char* OutList = NULL;
  bool AppendToOutput = false;
  
  bool Matrix3x3 = false;
  bool Transpose = false;
  bool Inverse = false;
  
  const char* InputFileName = NULL;
  
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Matrix to degrees of freedom" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Convert transformation matrix to degrees of freedom" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "mat2dof [options] < matrix" );

    typedef cmtk::CommandLine::Key Key;    
    cl.AddSwitch( Key( '3', "matrix3x3" ), &Matrix3x3, true, "Only input upper left 3x3 submatrix." );
    cl.AddSwitch( Key( 't', "transpose" ), &Transpose, true, "Transpose input matrix." );
    cl.AddSwitch( Key( 'i', "inverse" ), &Inverse, true, "Output inverse transformation." );
    cl.AddOption( Key( 'c', "center" ), &CenterStr, "Set center x,y,z for rotation and scaling." );
    cl.AddOption( Key( 'o', "offset" ), &OffsetStr, "Set offset dx,dy,dz for translation." );
    cl.AddOption( Key( 'x', "xlate" ), &TranslateStr, "Translate result relative by given vector." );
    cl.AddOption( Key( 'p', "pixel-size" ), &PixelSizeStr, "For matrices in pixel space, set pixel size." );

    cl.AddOption( Key( 'l', "list" ), &OutList, "Write output in list format to this archive" );
    cl.AddOption( Key( 'A', "append" ), &OutList, "Append output to this archive", &AppendToOutput );

    cl.Parse( argc, argv );

    InputFileName = cl.GetNextOptional();
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }

  cmtk::Types::Coordinate matrix[4][4];
  if ( InputFileName )
    {
    std::ifstream stream( InputFileName );
    if ( stream.good() )
      {
      ReadMatrix( matrix, stream, Matrix3x3 );
      }
    else
      {
      cmtk::StdErr << "ERROR: cannot open input file " << InputFileName << "\n";
      }
    }
  else
    {
    ReadMatrix( matrix, std::cin, Matrix3x3 );
    }
  
  if ( Transpose )
    {
    for ( int row = 0; row < 3; ++row )
      for ( int col = 0; col < row; ++col )
	{
	const cmtk::Types::Coordinate tmp = matrix[col][row];
	matrix[col][row] = matrix[row][col];
	matrix[row][col] = tmp;
	}
    }

  if ( PixelSizeStr )
    {
    float pixel[4] = { 1.0, 1.0, 1.0, 1.0 };
    if ( 3 == sscanf( PixelSizeStr, "%f,%f,%f", pixel+0, pixel+1, pixel+2 ) )
      {
      for ( int row = 0; row < 4; ++row )
	for ( int col = 0; col < 4; ++col )
	  matrix[col][row] /= pixel[row];
      }
    }

  cmtk::AffineXform::SmartPtr xform( new cmtk::AffineXform( matrix ) );

  if ( CenterStr )
    {
    float center[3];
    if ( 3 == sscanf( CenterStr, "%f,%f,%f", center+0, center+1, center+2 ) )
      {
      xform->ChangeCenter( cmtk::FixedVector<3,cmtk::Types::Coordinate>::FromPointer( center ) );
      }
    }
  
  if ( OffsetStr )
    {
    float offset[3];
    if ( 3 == sscanf( OffsetStr, "%f,%f,%f", offset+0, offset+1, offset+2 ) )
      {
      xform->Translate( xform->RotateScaleShear( cmtk::FixedVector<3,cmtk::Types::Coordinate>::FromPointer( offset ) ) );
      }
    }

  if ( TranslateStr )
    {
    float xlate[3];
    if ( 3 == sscanf( TranslateStr, "%f,%f,%f", xlate+0, xlate+1, xlate+2 ) )
      {
      xform->Translate( cmtk::FixedVector<3,cmtk::Types::Coordinate>::FromPointer( xlate ) );
      }
    }
  
  if ( OutList )
    {
    if ( AppendToOutput )
      {
      cmtk::ClassStreamOutput outStream( OutList, cmtk::ClassStreamOutput::MODE_APPEND );
      if ( Inverse )
	outStream << *xform->GetInverse();
      else
	outStream << *xform;
      }
    else
      {
      cmtk::StudyList studyList;
      
      studyList.AddStudy( Reference );
      studyList.AddStudy( Floating );

      if ( Inverse )
	{
	cmtk::AffineXform::SmartPtr inverse( xform->GetInverse() );
	studyList.AddXform( Reference, Floating, inverse );
	}
      else
	studyList.AddXform( Reference, Floating, xform );
      
      cmtk::ClassStreamStudyList::Write( OutList, &studyList );
      }
    }
  else
    {
    cmtk::CoordinateVector v;
    xform->GetParamVector( v );
    for ( unsigned int idx = 0; idx < v.Dim; ++idx )
      std::cout << "#" << idx << "\t" << v.Elements[idx] << "\n";
    }

  return 0;
}

#include "cmtkSafeMain"
