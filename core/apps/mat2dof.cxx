/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2014 SRI International
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
    for ( int j = 0; j < 3; ++j )
      {
      for ( int i = 0; i < 3; ++i )
	{
	stream >> matrix[j][i];
	}
      }
    matrix[3][3] = 1.0;
    }
  else
    {
    for ( int j = 0; j < 4; ++j )
      {
      for ( int i = 0; i < 4; ++i )
	{
	stream >> matrix[j][i];
	}
      }
    }
}

int
doMain( const int argc, const char* argv[] )
{
  const char* CenterStr = NULL;
  const char* OffsetStr = NULL;
  const char* TranslateStr = NULL;
  const char* PixelSizeStr = NULL;
  
  const char* OutList = NULL;
  const char* OutXform = NULL;
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
    cl.BeginGroup( "Input", "Input Options" );
    cl.AddSwitch( Key( '3', "matrix3x3" ), &Matrix3x3, true, "Only input upper left 3x3 submatrix." );
    cl.AddSwitch( Key( 't', "transpose" ), &Transpose, true, "Transpose input matrix." );
    cl.AddOption( Key( "pixel-size" ), &PixelSizeStr, "For matrices in pixel space (i.e., mapping from grid index to physical coordinate), set image grid pixel size." );
    cl.EndGroup();

    cl.BeginGroup( "Modifiers", "Transformation Modifiers" );
    cl.AddOption( Key( "center" ), &CenterStr, "Set center x,y,z for rotation and scaling." );
    cl.AddOption( Key( "offset" ), &OffsetStr, "Set offset dx,dy,dz for translation." );
    cl.AddOption( Key( "xlate" ), &TranslateStr, "Translate result relative by given vector." );
    cl.AddSwitch( Key( "inverse" ), &Inverse, true, "Output inverse transformation. Inversion happens after all other modifiers have been applied." );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( 'o', "output" ), &OutXform, "Write output transformation to this file" );
    cl.AddOption( Key( "list" ), &OutList, "Write output in list format to this archive (directory including 'registration' and 'studylist' files)." );
    cl.AddSwitch( Key( "append" ), &AppendToOutput, true, "Append to output file rather than create a new one and overwrite exisiting files." );
    cl.EndGroup();

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
    for ( int row = 0; row < 4; ++row )
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
    if ( 3 == sscanf( PixelSizeStr, "%15f,%15f,%15f", pixel+0, pixel+1, pixel+2 ) )
      {
      for ( int row = 0; row < 4; ++row )
	for ( int col = 0; col < 4; ++col )
	  matrix[col][row] /= pixel[row];
      }
    }

  // Create transformation from imported matrix
  cmtk::AffineXform::SmartPtr xform;
  try
    {
    xform = cmtk::AffineXform::SmartPtr( new cmtk::AffineXform( matrix ) );
    }
  catch ( const cmtk::AffineXform::MatrixType::SingularMatrixException& )
    {
    cmtk::StdErr << "ERROR: singular matrix encountered in cmtk::AffineXform constructor.\n";
    throw cmtk::ExitException( 1 );
    }

  if ( CenterStr )
    {
    float center[3];
    if ( 3 == sscanf( CenterStr, "%15f,%15f,%15f", center+0, center+1, center+2 ) )
      {
      xform->ChangeCenter( cmtk::FixedVector<3,cmtk::Types::Coordinate>::FromPointer( center ) );
      }
    }
  
  if ( OffsetStr )
    {
    float offset[3];
    if ( 3 == sscanf( OffsetStr, "%15f,%15f,%15f", offset+0, offset+1, offset+2 ) )
      {
      xform->Translate( xform->RotateScaleShear( cmtk::FixedVector<3,cmtk::Types::Coordinate>::FromPointer( offset ) ) );
      }
    }

  if ( TranslateStr )
    {
    float xlate[3];
    if ( 3 == sscanf( TranslateStr, "%15f,%15f,%15f", xlate+0, xlate+1, xlate+2 ) )
      {
      xform->Translate( cmtk::FixedVector<3,cmtk::Types::Coordinate>::FromPointer( xlate ) );
      }
    }

  // Invert transformation
  if ( Inverse )
    {
    try
      {
      xform = xform->GetInverse();
      }
    catch ( const cmtk::AffineXform::MatrixType::SingularMatrixException& )
      {
      cmtk::StdErr << "ERROR: affine transformation with singular matrix cannot be inverted\n";
      throw cmtk::ExitException( 1 );
      }
    }

  // Are we writing to a list archive?
  if ( OutList )
    {
    if ( AppendToOutput )
      {
      cmtk::ClassStreamOutput outStream( OutList, "registration", cmtk::ClassStreamOutput::MODE_APPEND );
      outStream << *xform;
      }
    else
      {
      cmtk::StudyList studyList;
      
      const char* strReference = "reference";
      const char* strFloating = "floating";
  
      studyList.AddStudy( strReference );
      studyList.AddStudy( strFloating );

      studyList.AddXform( strReference, strFloating, xform );
      
      cmtk::ClassStreamStudyList::Write( OutList, &studyList );
      }
    }

  // Are we writing to a single file?
  if ( OutXform )
    {
    cmtk::ClassStreamOutput outStream( OutXform, AppendToOutput ? cmtk::ClassStreamOutput::MODE_APPEND : cmtk::ClassStreamOutput::MODE_WRITE );
    outStream << *xform;
    }

  // No output files given: just dump parameters to console
  if ( !OutList && !OutXform )
    {
    cmtk::CoordinateVector v;
    if ( Inverse )
      {
      try 
	{
	xform->GetInverse()->GetParamVector( v );
	}
      catch ( const cmtk::AffineXform::MatrixType::SingularMatrixException& )
	{
	cmtk::StdErr << "ERROR: affine transformation with singular matrix cannot be inverted\n";
	throw cmtk::ExitException( 1 );
	}
      }
    else
      {
      xform->GetParamVector( v );
      }
    for ( unsigned int idx = 0; idx < v.Dim; ++idx )
      std::cout << "#" << idx << "\t" << v.Elements[idx] << "\n";
    }

  return 0;
}

#include "cmtkSafeMain"
