/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include "cmtkTypedStreamOutput.h"

#include <System/cmtkFileUtils.h>
#include <System/cmtkConsole.h>

#include <string.h>
#include <stdlib.h>
#include <limits.h>

#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

namespace
cmtk
{

/** \addtogroup IO */
//@{

TypedStreamOutput::TypedStreamOutput
( const std::string& filename, const Self::Mode mode )
{
  this->Open( filename, mode );
}

TypedStreamOutput::TypedStreamOutput
( const std::string& dir, const std::string& archive, const Self::Mode mode )
{
  this->Open( dir, archive, mode );
}

TypedStreamOutput
::~TypedStreamOutput()
{
  this->Close();
}

void 
TypedStreamOutput
::Open
( const std::string& dir, const std::string& archive, const Self::Mode mode )
{
  static char fname[PATH_MAX];
  
  // If "dir" parameter is empty, use current directory instead.
  if ( dir != "" ) 
    {
    if ( static_cast<size_t>( snprintf( fname, sizeof( fname ), "%s%c%s", dir.c_str(), CMTK_PATH_SEPARATOR, archive.c_str() ) ) >= sizeof( fname ) )
      {
      StdErr << "WARNING: length of path exceeds system PATH_MAX in TypedStreamOutput::Open and will be truncated.\n";
      }
    } 
  else 
    {
    if ( static_cast<size_t>( snprintf( fname, sizeof( fname ), "%s", archive.c_str() ) ) >= sizeof( fname ) )
      {
      StdErr << "WARNING: length of path exceeds system PATH_MAX in TypedStreamOutput::Open and will be truncated.\n";
      }
    }
  
#ifndef _MSC_VER
  // check if this is an existing directory; if it is, set access/modification time.
  // this is useful for dependency checking using file system timestamps.
  struct stat buf;
  if ( !stat( dir.c_str(), &buf ) && S_ISDIR( buf.st_mode ) )
    {
    utimes( dir.c_str(), NULL );
    }
#endif

  this->Open( fname, mode );
}

void 
TypedStreamOutput
::Open
( const std::string& filename, const Self::Mode mode )
{
  this->m_Status = Self::ERROR_NONE;

  this->Close();
  
  if ( mode != Self::MODE_WRITE && mode != Self::MODE_WRITE_ZLIB && mode != Self::MODE_APPEND ) 
    {
    this->m_Status = Self::ERROR_ARG;
    return;
    }
  
  if ( mode == Self::MODE_WRITE || mode == Self::MODE_WRITE_ZLIB ) 
    {
    if ( FileUtils::RecursiveMkPrefixDir( filename ) ) 
      {
      StdErr << "ERROR: could not recursively create path for \"" << filename << "\"\n";
      this->m_Status = Self::ERROR_SYSTEM;
      return;
      }
    }
  
  const char *modestr = "";
  switch ( mode ) 
    {
#ifdef _MSC_VER
    case Self::MODE_WRITE:	modestr = "wb"; break;
    case Self::MODE_WRITE_ZLIB:	modestr = "wb"; break;
    case Self::MODE_APPEND:	modestr = "ab"; break;
#else
    case Self::MODE_WRITE:	modestr = "w"; break;
    case Self::MODE_WRITE_ZLIB:	modestr = "w"; break;
    case Self::MODE_APPEND:	modestr = "a"; break;
#endif
    default: modestr = ""; break;
    }
  
  if ( mode == Self::MODE_WRITE_ZLIB )
    {
    const std::string gzName = filename + ".gz";
    GzFile = gzopen( gzName.c_str(), modestr );
    if ( ! GzFile ) 
      {
      StdErr << "ERROR: could not open gz file \"" << gzName << "\" with mode \"" << modestr << "\"\n";
      this->m_Status = Self::ERROR_SYSTEM;
      return;
      }
    }
  else
    {
    File = fopen( filename.c_str(), modestr );
    if ( ! File ) 
      {
      StdErr << "ERROR: could not open file \"" << filename << "\" with mode \"" << modestr << "\"\n";
      this->m_Status = Self::ERROR_SYSTEM;
      return;
      }
    }
  
  this->m_Mode = mode;
  switch ( this->m_Mode ) 
    {
    default: break;
    case Self::MODE_WRITE:
    case Self::MODE_WRITE_ZLIB:
      if ( GzFile )
	gzprintf( GzFile, "%s\n", Self::GetTypedStreamIdent() );
      else
	fprintf( File, "%s\n", Self::GetTypedStreamIdent() );
      break;
      
    case Self::MODE_APPEND:
      if ( GzFile ) 
	{
	if ( 0 == gztell( this->GzFile) )
	  gzprintf( GzFile, "%s\n", Self::GetTypedStreamIdent() );
	} 
      else
	if ( 0 == ftell( File) )
	  fprintf( File, "%s\n", Self::GetTypedStreamIdent() );
      break;
    }
}

void
TypedStreamOutput
::Close()
{
  if ( File || GzFile )
    {
    while ( ! LevelStack.empty() ) 
      {
      LevelStack.pop();
      int streamLevel = LevelStack.size();
      if ( GzFile ) 
	{
	for ( int level = 0; level < streamLevel; level++)
	  gzputs( GzFile, "\t" );
	gzputs( GzFile, "}\n" );
	} 
      else
	{
	for ( int level = 0; level < streamLevel; level++)
	  fputs( "\t", File );
	fputs( "}\n", File );
	}
      }
    } 

  if ( this->GzFile )
    {
    gzclose( this->GzFile );
    this->GzFile = NULL;
    }
  
  if ( this->File )
    {
    fclose( this->File );
    this->File = NULL;
    }
  
  this->m_Status = Self::ERROR_NONE;
  SplitPosition = NULL;
}

TypedStreamOutput::Condition 
TypedStreamOutput
::Begin
( const std::string& section )
{
  if ( !File && !GzFile)
    {
    this->m_Status = Self::ERROR_INVALID;
    return Self::CONDITION_ERROR;
    }

  int streamLevel = LevelStack.size();
  if ( GzFile ) 
    {
    for ( int level = 0; level < streamLevel; level++)
      gzputs( GzFile, "\t" );
    
    gzprintf( GzFile, "%s {\n", section.c_str() );
    }
  else
    {
    for ( int level = 0; level < streamLevel; level++)
      fputs( "\t", File );
    
    fprintf( File, "%s {\n", section.c_str() );
    }
  
  if ( GzFile )
    LevelStack.push( gztell( GzFile ) );
  else
    LevelStack.push( ftell( File ) );
  
  return Self::CONDITION_OK;
}

TypedStreamOutput::Condition
TypedStreamOutput
::End
( const bool flush )
{
  if ( ! File && ! GzFile )
    {
    this->m_Status = Self::ERROR_INVALID;
    return Self::CONDITION_ERROR;
    }
  
  int streamLevel = LevelStack.size();
  if ( streamLevel == 0 ) 
    {
    // end without begin
    this->m_Status = Self::ERROR_LEVEL;
    return Self::CONDITION_ERROR;
    }
  
  LevelStack.pop();
  if ( GzFile ) 
    {
    for ( int level = 0; level < streamLevel-1; level++ )
      gzputs( GzFile, "\t" );
    gzputs( GzFile, "}\n" );
    }
  else 
    {
    for ( int level = 0; level < streamLevel-1; level++ )
      fputs( "\t", File );
    fputs( "}\n", File ); 
    }
  
  if ( flush ) 
    {
    fflush( File );
    }
  
  return Self::CONDITION_OK;
}
  
TypedStreamOutput::Condition 
TypedStreamOutput
::WriteBool( const char* key, const bool value )
{
  int currentLevel = LevelStack.size();
  if ( GzFile ) 
    {
    for ( int level = 0; level < currentLevel; level++ )
      gzputs( GzFile, "\t" );
    gzprintf( GzFile, "%s %s\n", key, (value) ? "yes" : "no");
    } 
  else
    { 
    for ( int level = 0; level < currentLevel; level++ )
      fputs( "\t", File );
    fprintf( File, "%s %s\n", key, (value) ? "yes" : "no");
    }
  return Self::CONDITION_OK;
}

TypedStreamOutput::Condition
TypedStreamOutput
::WriteInt( const char* key, const int value )
{
  int currentLevel = LevelStack.size();
  if ( GzFile ) 
    {
    for ( int level = 0; level < currentLevel; level++ )
      gzputs( GzFile, "\t" );
    gzprintf( GzFile, "%s %d\n", key, value );
    } 
  else
    { 
    for ( int level = 0; level < currentLevel; level++ )
      fputs( "\t", File );
    fprintf( File, "%s %d\n", key, value );
    }
  return Self::CONDITION_OK;
}

TypedStreamOutput::Condition
TypedStreamOutput
::WriteFloat( const char* key, const float value )
{
  int currentLevel = LevelStack.size();
  if ( GzFile ) 
    {
    for ( int level = 0; level < currentLevel; level++ )
      gzputs( GzFile, "\t" );
    gzprintf( GzFile, "%s %.*f\n", key, PrecisionFloat, value );
    } 
  else
    {
    for ( int level = 0; level < currentLevel; level++ )
      fputs( "\t", File );
    fprintf( File, "%s %.*f\n", key, PrecisionFloat, value );
    }
  return Self::CONDITION_OK;
}

TypedStreamOutput::Condition
TypedStreamOutput
::WriteDouble( const char* key, const double value )
{
  int currentLevel = LevelStack.size();
  if ( GzFile ) 
    {
    for ( int level = 0; level < currentLevel; level++ )
      gzputs( GzFile, "\t" );
    gzprintf( this->GzFile, "%s %.*f\n", key, PrecisionDouble, value );
    } 
  else
    {
    for ( int level = 0; level < currentLevel; level++ )
      fputs( "\t", File );
    fprintf( File, "%s %.*f\n", key, PrecisionDouble, value );
    }
  return Self::CONDITION_OK;
}

TypedStreamOutput::Condition
TypedStreamOutput
::WriteString( const char* key, const std::string& value )
{
  return this->WriteString( key, value.c_str() );
}

TypedStreamOutput::Condition
TypedStreamOutput
::WriteString( const char* key, const char* value )
{
  char *buffer = Buffer;
  
  const char *strValue = (value) ? value : "";
  while (*strValue) 
    {
    if (*strValue == '\\') 
      {
      *buffer++ = '\\';
      *buffer++ = *strValue++;
      continue;
      }
    if (*strValue == '\"') 
      {
      *buffer++ = '\\';
      *buffer++ = *strValue++;
      continue;
      }
    if (*strValue == '\n') 
      {
      *buffer++ = '\\';
      *buffer++ = 'n';
      strValue++;
      continue;
      }
    *buffer++ = *strValue++;
    }
  
  *buffer++ = '\0';
  
  int currentLevel = LevelStack.size();
  if ( GzFile ) 
    {
    for ( int level = 0; level < currentLevel; level++)
      gzputs( GzFile, "\t" );
    gzprintf( GzFile, "%s \"%s\"\n", key, Buffer);
    } 
  else
    {
    for ( int level = 0; level < currentLevel; level++)
      fputs( "\t", File );
    fprintf( File, "%s \"%s\"\n", key, Buffer);
    }
  
  return Self::CONDITION_OK;
}

TypedStreamOutput::Condition
TypedStreamOutput
::WriteIntArray( const char* key, const int* array, const int size, const int valuesPerLine )
{
  if ( !array || size < 1) 
    {
    this->m_Status = Self::ERROR_ARG;
    return Self::CONDITION_ERROR;
    }
  
  int currentLevel = LevelStack.size();
  if ( GzFile ) 
    {
    for ( int level = 0; level < currentLevel; level++)
      gzputs( GzFile, "\t" );
    
    gzprintf( this->GzFile, "%s ", key );
    
    for ( int i = 0; i < size; i++) 
      {
      if (i && (i % valuesPerLine) == 0) 
	{
	gzprintf( GzFile, "\n\t");
	for ( int level = 0; level < currentLevel; level++ )
	  gzputs( GzFile, "\t" );
	}
      gzprintf( GzFile, "%d ", array[i] );
      }
    
    gzputs( GzFile, "\n" );
    } 
  else
    {
    for ( int level = 0; level < currentLevel; level++)
      fputs( "\t", File );
    
    fprintf( File, "%s ", key );
    
    for ( int i = 0; i < size; i++) 
      {
      if (i && (i % valuesPerLine) == 0) 
	{
	fprintf( File, "\n\t");
	for ( int level = 0; level < currentLevel; level++ )
	  fputs( "\t", File );
	}
      fprintf( File, "%d ", array[i] );
      }
    
    fputs( "\n", File);
    }
  
  return Self::CONDITION_OK;
}

TypedStreamOutput::Condition
TypedStreamOutput
::WriteBoolArray( const char* key, const byte* array, const int size, const int valuesPerLine )
{
  if ( !array || size < 1) 
    {
    this->m_Status = Self::ERROR_ARG;
    return Self::CONDITION_ERROR;
    }
  
  int currentLevel = LevelStack.size();
  if ( GzFile ) 
    {
    for ( int level = 0; level < currentLevel; level++)
      gzputs( GzFile, "\t" );
    
    gzprintf( GzFile, "%s ", key );
    
    for ( int i = 0; i < size; i++) 
      {
      if (i && (i % valuesPerLine) == 0) 
	{
	gzprintf( GzFile, "\n\t");
	for ( int level = 0; level < currentLevel; level++ )
	  gzputs( GzFile, "\t" );
	}
      gzprintf( GzFile, "%d", (array[i/8]>>(i%8))&1 );
      }
    
    gzputs( GzFile, "\n" );
    }
  else
    {
    for ( int level = 0; level < currentLevel; level++)
      fputs( "\t", File );
    
    fprintf( File, "%s ", key );
    
    for ( int i = 0; i < size; i++) 
      {
      if (i && (i % valuesPerLine) == 0) 
	{
	fprintf( File, "\n\t");
	for ( int level = 0; level < currentLevel; level++ )
	  fputs( "\t", File );
	}
      fprintf( File, "%d", (array[i/8]>>(i%8))&1 );
      }
    
    fputs( "\n", File );
    }
  
  return Self::CONDITION_OK;
}

TypedStreamOutput::Condition
TypedStreamOutput
::WriteFloatArray( const char* key, const float* array, const int size, const int valuesPerLine )
{
  if ( !array || size < 1) 
    {
    this->m_Status = Self::ERROR_ARG;
    return Self::CONDITION_ERROR;
    }
  
  int currentLevel = LevelStack.size();
  if ( GzFile ) 
    {
    for ( int level = 0; level < currentLevel; level++)
      gzputs( GzFile, "\t" );
    
    gzprintf( GzFile, "%s ", key );
    
    for ( int i = 0; i < size; i++ ) 
      {
      if (i && (i % valuesPerLine) == 0) 
	{
	gzprintf( GzFile, "\n\t");
	for ( int level = 0; level < currentLevel; level++ )
	  gzputs( GzFile, "\t" );
	}
      gzprintf( GzFile, "%.*g ", PrecisionFloat, array[i] );
      }
    
    gzprintf( GzFile, "\n" );
    } 
  else
    {
    for ( int level = 0; level < currentLevel; level++)
      fputs( "\t", File );
    
    fprintf( File, "%s ", key );
    
    for ( int i = 0; i < size; i++ ) 
      {
      if (i && (i % valuesPerLine) == 0) 
	{
	fprintf( File, "\n\t");
	for ( int level = 0; level < currentLevel; level++ )
	  fputs( "\t", File );
	}
      fprintf( File, "%.*g ", PrecisionFloat, array[i] );
      }
    
    fprintf( File, "\n" );
    }
  
  return Self::CONDITION_OK;
}

TypedStreamOutput::Condition
TypedStreamOutput
::WriteDoubleArray( const char* key, const double* array, const int size, const int valuesPerLine )
{
  if ( !array || size < 1) 
    {
    this->m_Status = Self::ERROR_ARG;
    return Self::CONDITION_ERROR;
    }
  
  int currentLevel = LevelStack.size();
  if ( GzFile )
    {
    for ( int level = 0; level < currentLevel; level++)
      gzputs( GzFile, "\t" );
    
    gzprintf( GzFile, "%s ", key );
    
    for ( int i = 0; i < size; i++ ) 
      {
      if (i && (i % valuesPerLine) == 0) 
	{
	gzprintf( GzFile, "\n\t");
	for ( int level = 0; level < currentLevel; level++ )
	  gzputs( GzFile, "\t" );
	}
      gzprintf( GzFile, "%.*g ", PrecisionDouble, array[i] );
      }
    
    gzprintf( GzFile, "\n" );
    } 
  else
    {
    for ( int level = 0; level < currentLevel; level++)
      fputs( "\t", File );
    
    fprintf( File, "%s ", key );
    
    for ( int i = 0; i < size; i++ ) 
      {
      if (i && (i % valuesPerLine) == 0) 
	{
	fprintf( File, "\n\t");
	for ( int level = 0; level < currentLevel; level++ )
	  fputs( "\t", File );
	}
      fprintf( File, "%.*g ", PrecisionDouble, array[i] );
      }
    
    fprintf( File, "\n" );
    }
  
  return Self::CONDITION_OK;
}

} // namespace cmtk
