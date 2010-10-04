/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkTypedStream.h"

#include <System/cmtkFileUtils.h>
#include <System/cmtkConsole.h>

#include <string.h>
#include <stdlib.h>
#include <cstdarg>
#include <limits.h>

#ifdef HAVE_MALLOC_H
#  include <malloc.h>
#endif

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

TypedStream
::TypedStream()
{
  this->InitInternals();
}

TypedStream
::TypedStream
( const char *filename, const TypedStreamMode mode )
{
  this->InitInternals();
  this->Open( filename, mode );
}

TypedStream
::TypedStream
( const char *dir, const char *archive, const TypedStreamMode mode )
{
  this->InitInternals();
  this->Open( dir, archive, mode );
}

void
TypedStream
::InitInternals()
{
  File = NULL;
  GzFile = NULL;

  PrecisionFloat = 6;
  PrecisionDouble = 10;

  Status = TYPEDSTREAM_ERROR_NONE;
  DebugFlag = TYPEDSTREAM_DEBUG_OFF;

  SplitPosition = NULL;
}

TypedStream
::~TypedStream()
{
  this->Close();
  this->InitInternals();
}

void 
TypedStream
::Open
( const char *dir, const char *archive, const TypedStreamMode mode )
{
  static char fname[PATH_MAX];
  
  // If "dir" parameter is NULL or empty, use current directory instead.
  if ( dir && *dir ) 
    {
    if ( static_cast<size_t>( snprintf( fname, sizeof( fname ), "%s/%s", dir, archive ) ) >= sizeof( fname ) )
      {
      StdErr << "WARNING: length of path exceeds system PATH_MAX in TypedStream::Open and will be truncated.\n";
      }
    } 
  else 
    {
    if ( static_cast<size_t>( snprintf( fname, sizeof( fname ), "%s", archive ) ) >= sizeof( fname ) )
      {
      StdErr << "WARNING: length of path exceeds system PATH_MAX in TypedStream::Open and will be truncated.\n";
      }
    }
  
  if ( mode == TYPEDSTREAM_WRITE || mode == TYPEDSTREAM_WRITE_ZLIB || mode == TYPEDSTREAM_APPEND ) 
    {
#ifndef _MSC_VER
    // check if this is an existing directory; if it is, set access/modification time.
    // this is useful for dependency checking using file system timestamps.
    struct stat buf;
    if ( !stat( dir, &buf ) && S_ISDIR( buf.st_mode ) )
      {
      utimes( dir, NULL );
      }
#endif
    }    

  this->Open( fname, mode );
}

void 
TypedStream
::Open
( const char *filename, const TypedStreamMode mode )
{
  Status = TYPEDSTREAM_ERROR_NONE;

  this->Close();
  
  if ( !filename) 
    {
    Status = TYPEDSTREAM_ERROR_ARG;
    return;
    }
  
  if ( mode != TYPEDSTREAM_READ && mode != TYPEDSTREAM_WRITE && mode != TYPEDSTREAM_WRITE_ZLIB && mode != TYPEDSTREAM_APPEND ) 
    {
    Status = TYPEDSTREAM_ERROR_ARG;
    return;
    }
  
  if ( mode == TYPEDSTREAM_WRITE || mode == TYPEDSTREAM_WRITE_ZLIB ) 
    {
    if ( FileUtils::RecursiveMkPrefixDir( filename ) ) 
      {
      Status = TYPEDSTREAM_ERROR_SYSTEM;
      return;
      }
    }
  
  const char *modestr = "";
  switch ( mode ) 
    {
#ifdef _MSC_VER
    case TYPEDSTREAM_READ:	modestr = "rb"; break;
    case TYPEDSTREAM_WRITE:	modestr = "wb"; break;
    case TYPEDSTREAM_WRITE_ZLIB:	modestr = "wb"; break;
    case TYPEDSTREAM_APPEND:	modestr = "ab"; break;
#else
    case TYPEDSTREAM_READ:	modestr = "r"; break;
    case TYPEDSTREAM_WRITE:	modestr = "w"; break;
    case TYPEDSTREAM_WRITE_ZLIB:	modestr = "w"; break;
    case TYPEDSTREAM_APPEND:	modestr = "a"; break;
#endif
    }
  
  if ( ! ( File = fopen( filename, modestr ) ) ) 
    {
    char gzName[PATH_MAX];
    snprintf( gzName, sizeof( gzName ), "%s.gz", filename );
    GzFile = gzopen( gzName, modestr );
    if ( ! GzFile ) 
      {
      Status = TYPEDSTREAM_ERROR_SYSTEM;
      return;
      }
    }
  
  Mode = mode;
  switch ( Mode ) 
    {
    case TYPEDSTREAM_READ:
      if ( GzFile ) 
	{
	if (! gzgets( GzFile, Buffer, TYPEDSTREAM_LIMIT_BUFFER ) ) 
	  {
	  Status = TYPEDSTREAM_ERROR_FORMAT;
	  gzclose( GzFile );
	  return;
	  }
	} 
      else
	if (! fgets( Buffer, TYPEDSTREAM_LIMIT_BUFFER, File ) ) 
	  {
	  Status = TYPEDSTREAM_ERROR_FORMAT;
	  fclose( File );
	  return;
	  }
      
      if (Buffer[0] != '!' && Buffer[0] != '#') 
	{
	Status = TYPEDSTREAM_ERROR_FORMAT;
	if ( GzFile )
	  gzclose( GzFile );
	else
	  fclose( File );
	return;
	}
      int releaseMinor, releaseMajor;
      if (2 != sscanf( Buffer+1, " TYPEDSTREAM %d.%d", &releaseMajor, &releaseMinor)) 
	{
	Status = TYPEDSTREAM_ERROR_FORMAT;
	if ( GzFile )
	  gzclose( GzFile );
	else
	  fclose( File );
	return;
	}
      if ( releaseMajor == 1 && releaseMinor == 0) 
	{
	/* Release 1.0 */
	} 
      else
	if (releaseMajor == 1 && releaseMinor == 1) 
	  {
	  /* Release 1.1 */
	  } 
	else 
	  {
	  /* Unknown Release */
	  Status = TYPEDSTREAM_ERROR_FORMAT;
	  if ( GzFile )
	    gzclose( GzFile );
	  else
	    fclose( File );
	  return;
	  }
      break;
      
    case TYPEDSTREAM_WRITE:
    case TYPEDSTREAM_WRITE_ZLIB:
      if ( GzFile )
	gzprintf( GzFile, "%s\n", GetTypedStreamIdent() );
      else
	fprintf( File, "%s\n", GetTypedStreamIdent() );
      break;
      
    case TYPEDSTREAM_APPEND:
      if ( GzFile ) 
	{
	if ( 0 == gztell( File) )
	  gzprintf( GzFile, "%s\n", GetTypedStreamIdent() );
	} 
      else
	if ( 0 == ftell( File) )
	  fprintf( File, "%s\n", GetTypedStreamIdent() );
      break;
    }
}

void
TypedStream
::Close()
{
  if ( File || GzFile )
    {
    if ( Mode == TYPEDSTREAM_WRITE || Mode == TYPEDSTREAM_APPEND ) 
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
    else 
      {
      while ( ! LevelStack.empty() ) 
	{
	LevelStack.pop();      
	}
      }
    
    if ( GzFile ) 
      {
      gzclose( GzFile );
      GzFile = NULL;
      }
    if ( File ) 
      {
      fclose( File );
      File = NULL;
      }
    }
  
  Status = TYPEDSTREAM_ERROR_NONE;
  SplitPosition = NULL;
}

TypedStreamCondition 
TypedStream
::Begin
( const char* section )
{
  if ( !File && !GzFile)
    {
    Status = TYPEDSTREAM_ERROR_INVALID;
    return TYPEDSTREAM_ERROR;
    }

  if ( (Mode != TYPEDSTREAM_READ) && !section) 
    {
    Status = TYPEDSTREAM_ERROR_ARG;
    return TYPEDSTREAM_ERROR;
    }
  
  if ( Mode == TYPEDSTREAM_READ ) 
    {
    if ( GzFile )
      gzseek( GzFile, LevelStack.top(), SEEK_SET );
    else
      fseek( File, LevelStack.top(), SEEK_SET );
    return TYPEDSTREAM_OK;
    } 
  else 
    {
    int streamLevel = LevelStack.size();
    if ( GzFile ) 
      {
      for ( int level = 0; level < streamLevel; level++)
	gzputs( GzFile, "\t" );

      gzprintf( GzFile, "%s {\n", section );
      }
    else
      {
      for ( int level = 0; level < streamLevel; level++)
	fputs( "\t", File );
      
      fprintf( File, "%s {\n", section );
      }
    
    if ( GzFile )
      LevelStack.push( gztell( GzFile ) );
    else
      LevelStack.push( ftell( File ) );
    }
  
  return TYPEDSTREAM_OK;
}

TypedStreamCondition
TypedStream
::End
( const bool flush )
{
  if ( ! File && ! GzFile )
    {
    Status = TYPEDSTREAM_ERROR_INVALID;
    return TYPEDSTREAM_ERROR;
    }
  
  if ( Mode == TYPEDSTREAM_READ ) 
    {
    if ( LevelStack.size() == 0 ) 
      {
      // end without begin
      Status = TYPEDSTREAM_ERROR_LEVEL;
      return TYPEDSTREAM_ERROR;
      }

    TypedStreamToken token;
    int currentLevel = 1;
    while ( currentLevel && (TYPEDSTREAM_EOF != ( token = this->ReadLineToken() ) ) ) 
      {
      if ( token == TYPEDSTREAM_BEGIN ) 
	{
	this->DebugOutput( "Skipping section %s at level %d.", BufferKey, currentLevel );
	++currentLevel;
	}
      else
	if ( token == TYPEDSTREAM_END ) 
	  {
	  this->DebugOutput( "Leaving section %d.", currentLevel );
	  --currentLevel;
	  }
      }
    
    LevelStack.pop();
    } 
  else 
    {
    int streamLevel = LevelStack.size();
    if ( streamLevel == 0 ) 
      {
      // end without begin
      Status = TYPEDSTREAM_ERROR_LEVEL;
      return TYPEDSTREAM_ERROR;
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
    }
  
  if ( flush ) 
    {
    fflush( File );
    }
  
  return TYPEDSTREAM_OK;
}
  
TypedStreamCondition
TypedStream
::Seek
( const char* section, const bool forward )
{
  if ( ! File && ! GzFile )
    {
    Status = TYPEDSTREAM_ERROR_INVALID;
    return TYPEDSTREAM_ERROR;
    }

  if ( ! section ) 
    {
    Status = TYPEDSTREAM_ERROR_ARG;
    return TYPEDSTREAM_ERROR;
    }
  if ( Mode != TYPEDSTREAM_READ ) 
    {
    Status = TYPEDSTREAM_ERROR_MODE;
    return TYPEDSTREAM_ERROR;
    }
  
  unsigned initialLevel = LevelStack.size();

  if ( ! forward )
    {
    if ( GzFile ) 
      {
      if ( initialLevel ) 
	{
	gzseek( GzFile, LevelStack.top(), SEEK_SET );
	} 
      else 
	{
	gzseek( GzFile, 0, SEEK_SET);
	}
      } 
    else
      if ( initialLevel ) 
	{
	fseek( File, LevelStack.top(), SEEK_SET );
	} 
      else 
	{
	fseek( File, 0, SEEK_SET);
	}
    }

  unsigned currentLevel = initialLevel;
  
  this->DebugOutput( "Seeking section %s from level %d.", section, initialLevel );
  TypedStreamToken token;
  while ( TYPEDSTREAM_EOF != ( token = this->ReadLineToken() ) ) 
    {
    if ( token == TYPEDSTREAM_BEGIN ) 
      {
      this->DebugOutput( "Enter section %s at level %d.", BufferKey, currentLevel );
      if ( this->StringCmp( BufferKey, section ) == 0 ) 
	{
	if ( currentLevel == LevelStack.size() ) 
	  {
	  if ( GzFile )
	    LevelStack.push( gztell( GzFile ) );
	  else
	    LevelStack.push( ftell( File ) );
	  return TYPEDSTREAM_OK;
	  }
	if ( currentLevel == LevelStack.size()-1 ) 
	  {
	  LevelStack.pop();
	  if ( GzFile )
	    LevelStack.push( gztell( GzFile ) );
	  else
	    LevelStack.push( ftell( File ) );
	  return TYPEDSTREAM_OK;
	  }
	}
      ++ currentLevel;
      }
    if ( token == TYPEDSTREAM_END ) 
      {
      this->DebugOutput( "Leaving section %d.", currentLevel );
      if ( ! currentLevel ) 
	{
	Status = TYPEDSTREAM_ERROR_LEVEL;
	return TYPEDSTREAM_ERROR;
	}
      if ( currentLevel < initialLevel ) 
	{
	Status = TYPEDSTREAM_ERROR_NONE;
	return TYPEDSTREAM_ERROR;
	}
      -- currentLevel;
      }
    }
  
  this->DebugOutput( "Section %s not found.", section );
  Status = TYPEDSTREAM_ERROR_NONE;
  return TYPEDSTREAM_ERROR;
}

TypedStreamCondition
TypedStream
::Rewind()
{
  if ( ! File && ! GzFile )
    {
    Status = TYPEDSTREAM_ERROR_INVALID;
    return TYPEDSTREAM_ERROR;
    }
  if ( Mode != TYPEDSTREAM_READ ) 
    {
    Status = TYPEDSTREAM_ERROR_MODE;
    return TYPEDSTREAM_ERROR;
    }
  
  if ( LevelStack.size() )
    LevelStack.pop();
  
  if ( LevelStack.size() == 0 ) 
    {
    if ( GzFile ) 
      {
      if ( -1 == gzseek( GzFile, 0, SEEK_SET ) ) 
	{
	Status = TYPEDSTREAM_ERROR_SYSTEM;
	return TYPEDSTREAM_ERROR;
	}
      } 
    else
      if ( -1 == fseek( File, 0, SEEK_SET ) )
	{
	Status = TYPEDSTREAM_ERROR_SYSTEM;
	return TYPEDSTREAM_ERROR;
	}
    } 
  else 
    {
    if ( GzFile ) 
      {
      if (-1 == gzseek( GzFile, LevelStack.top(), SEEK_SET)) 
	{
	Status = TYPEDSTREAM_ERROR_SYSTEM;
	return TYPEDSTREAM_ERROR;
	}
      }
    else
      if (-1 == fseek( File, LevelStack.top(), SEEK_SET)) 
	{
	Status = TYPEDSTREAM_ERROR_SYSTEM;
	return TYPEDSTREAM_ERROR;
	}
    }
  return TYPEDSTREAM_OK;
}

bool
TypedStream
::ReadBool( const char* key, const bool defaultValue, const bool forward )
{
  int value;
  
  if ( this->GenericReadArray( key, TYPEDSTREAM_TYPE_BOOL, &value, 1, forward ) != TYPEDSTREAM_OK )
    if ( this->GenericReadArray( key, TYPEDSTREAM_TYPE_INT, &value, 1, forward ) != TYPEDSTREAM_OK )
      return defaultValue;
  
  return (value != 0);
}

TypedStreamCondition 
TypedStream
::ReadBoolArray( const char* key, byte *const array, const int size, const bool forward )
{
  return this->GenericReadArray( key, TYPEDSTREAM_TYPE_BINARYBOOL, array, size, forward );
}

int 
TypedStream
::ReadInt( const char* key, const int defaultValue, const bool forward )
{
  int value = defaultValue;
  if ( this->GenericReadArray( key, TYPEDSTREAM_TYPE_INT, &value, 1, forward ) != TYPEDSTREAM_OK )
    return defaultValue;
  
  return value;
}

TypedStreamCondition 
TypedStream
::ReadIntArray( const char* key, int *const array, const int size, const bool forward )
{
  return this->GenericReadArray( key, TYPEDSTREAM_TYPE_INT, array, size, forward );
}

float
TypedStream
::ReadFloat( const char* key, const float defaultValue, const bool forward )
{
  float value = defaultValue;
  if ( this->GenericReadArray( key, TYPEDSTREAM_TYPE_FLOAT, &value, 1, forward ) != TYPEDSTREAM_OK )
    return defaultValue;
  
  return value;
}

TypedStreamCondition
TypedStream
::ReadFloatArray( const char* key, float *const array, const int size, const bool forward )
{
  return this->GenericReadArray( key, TYPEDSTREAM_TYPE_FLOAT, array, size, forward );
}

double
TypedStream
::ReadDouble( const char* key, const double defaultValue, const bool forward )
{
  double value = defaultValue;
  if ( this->GenericReadArray( key, TYPEDSTREAM_TYPE_DOUBLE, &value, 1, forward ) != TYPEDSTREAM_OK )
    return defaultValue;
  
  return value;
}

TypedStreamCondition
TypedStream
::ReadDoubleArray( const char* key, double *const array, const int size, const bool forward )
{
  return this->GenericReadArray( key, TYPEDSTREAM_TYPE_DOUBLE, array, size, forward );
}

char*
TypedStream
::ReadString( const char* key, const char *defaultValue, const bool forward )
{
  char *value;
  if ( this->GenericReadArray( key, TYPEDSTREAM_TYPE_STRING, &value, 1, forward ) != TYPEDSTREAM_OK )
    {
    if ( defaultValue )
      return strdup( defaultValue );
    else
      return NULL;
    }
  
  return value;
}

TypedStreamCondition 
TypedStream
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
  return TYPEDSTREAM_OK;
}

TypedStreamCondition
TypedStream
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
  return TYPEDSTREAM_OK;
}

TypedStreamCondition
TypedStream
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
  return TYPEDSTREAM_OK;
}

TypedStreamCondition
TypedStream
::WriteDouble( const char* key, const double value )
{
  int currentLevel = LevelStack.size();
  if ( GzFile ) 
    {
    for ( int level = 0; level < currentLevel; level++ )
      gzputs( GzFile, "\t" );
    gzprintf( File, "%s %.*f\n", key, PrecisionDouble, value );
    } 
  else
    {
    for ( int level = 0; level < currentLevel; level++ )
      fputs( "\t", File );
    fprintf( File, "%s %.*f\n", key, PrecisionDouble, value );
    }
  return TYPEDSTREAM_OK;
}

TypedStreamCondition
TypedStream
::WriteString( const char* key, const std::string& value )
{
  return this->WriteString( key, value.c_str() );
}

TypedStreamCondition
TypedStream
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
  
  return TYPEDSTREAM_OK;
}

TypedStreamCondition 
TypedStream
::WriteComment( const char* fmt, ... )
{
  if ( this->Mode != TYPEDSTREAM_WRITE ) 
    {
    Status = TYPEDSTREAM_ERROR_MODE;
    return TYPEDSTREAM_ERROR;
    }

  static char buffer[1024];

  va_list args;
  va_start(args, fmt);
  vsnprintf( buffer, sizeof( buffer ), fmt, args );
  va_end(args);

  if ( this->GzFile )
    gzprintf( GzFile, "! %s\n", buffer );
  else
    fprintf( File, "! %s\n", buffer );

  return TYPEDSTREAM_OK;
}

TypedStreamCondition 
TypedStream
::WriteComment( int argc, char* argv[] )
{
  return this->WriteComment( argc, const_cast<const char**>( argv ) );
}

TypedStreamCondition 
TypedStream
::WriteComment( const int argc, const char* argv[] )
{
  if ( this->Mode != TYPEDSTREAM_WRITE ) 
    {
    Status = TYPEDSTREAM_ERROR_MODE;
    return TYPEDSTREAM_ERROR;
    }

  if ( this->GzFile )
    gzprintf( GzFile, "! " );
  else
    fprintf( File, "! " );

  for ( int i = 0; i < argc; ++i )
    {
  if ( this->GzFile )
    gzprintf( GzFile, "%s ", argv[i] );
  else
    fprintf( File, "%s ", argv[i] );
    }

  if ( this->GzFile )
    gzputs( GzFile, "\n" );
  else
    fputs( "\n", File );

  return TYPEDSTREAM_OK;
}

TypedStreamCondition
TypedStream
::WriteIntArray( const char* key, const int* array, const int size, const int valuesPerLine )
{
  if ( !array || size < 1) 
    {
    Status = TYPEDSTREAM_ERROR_ARG;
    return TYPEDSTREAM_ERROR;
    }
  
  int currentLevel = LevelStack.size();
  if ( GzFile ) 
    {
    for ( int level = 0; level < currentLevel; level++)
      gzputs( GzFile, "\t" );
    
    gzprintf( File, "%s ", key );
    
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
  
  return TYPEDSTREAM_OK;
}

TypedStreamCondition
TypedStream
::WriteBoolArray( const char* key, const byte* array, const int size, const int valuesPerLine )
{
  if ( !array || size < 1) 
    {
    Status = TYPEDSTREAM_ERROR_ARG;
    return TYPEDSTREAM_ERROR;
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
  
  return TYPEDSTREAM_OK;
}

TypedStreamCondition
TypedStream
::WriteFloatArray( const char* key, const float* array, const int size, const int valuesPerLine )
{
  if ( !array || size < 1) 
    {
    Status = TYPEDSTREAM_ERROR_ARG;
    return TYPEDSTREAM_ERROR;
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
  
  return TYPEDSTREAM_OK;
}

TypedStreamCondition
TypedStream
::WriteDoubleArray( const char* key, const double* array, const int size, const int valuesPerLine )
{
  if ( !array || size < 1) 
    {
    Status = TYPEDSTREAM_ERROR_ARG;
    return TYPEDSTREAM_ERROR;
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
  
  return TYPEDSTREAM_OK;
}

TypedStreamCondition
TypedStream
::GenericReadArray( const char * key, const int type, void *const array, const int arraySize, const bool forward )
{
  if (!array || arraySize < 1) 
    {
    Status = TYPEDSTREAM_ERROR_ARG;
    return TYPEDSTREAM_ERROR;
    }
  
  unsigned currentLevel = LevelStack.size();
  if ( ! forward ) 
    {
    if ( GzFile ) 
      {
      if ( currentLevel )
	gzseek( GzFile, LevelStack.top(), SEEK_SET );
      else
	gzseek( GzFile, 0, SEEK_SET );
      }
    else
      if ( currentLevel )
	fseek( File, LevelStack.top(), SEEK_SET );
      else
	fseek( File, 0, SEEK_SET );
    }
  
  int line;
  char *buffer;
  while ( TYPEDSTREAM_EOF != ( line = this->ReadLineToken() ) ) 
    {
    if ( line == TYPEDSTREAM_KEY) 
      {
      if ( (currentLevel == LevelStack.size()) && (this->StringCmp( BufferKey, key )) == 0 ) 
	{
	int i = 0;
	switch (type) 
	  {
	  case TYPEDSTREAM_TYPE_INT: 
	  {
	  int *arrayInt = static_cast<int*>( array ); 
	  do
	    {
	    buffer = strtok( BufferValue, "\t\n " );
	    while ( buffer ) 
	      {
	      if ( *buffer == '\"' ) 
		{
		Status = TYPEDSTREAM_ERROR_TYPE;
		return TYPEDSTREAM_ERROR;
		}
	      if ( i >= arraySize )
		return TYPEDSTREAM_OK;
	      arrayInt[i++] = atoi(buffer);
	      buffer = strtok( NULL, "\t\n " );
	      }
	    } while ( i < arraySize && TYPEDSTREAM_VALUE == this->ReadLineToken() );
	  if ( i < arraySize) 
	    {
	    Status = TYPEDSTREAM_ERROR_ARG;
	    return TYPEDSTREAM_ERROR;
	    }
	  break;
	  }
	  case TYPEDSTREAM_TYPE_BOOL: 
	  {
	  int *arrayInt = static_cast<int*>( array ); 
	  if (arraySize != 1) 
	    {
	    Status = TYPEDSTREAM_ERROR_ARG;
	    return TYPEDSTREAM_ERROR;
	    }
	  do 
	    {
	    buffer = strtok( BufferValue, "\t\n " );
	    while ( buffer ) 
	      {
	      if ( i >= arraySize )
		return TYPEDSTREAM_OK;
	      if ((buffer[0] == 'y' || buffer[0] == 'Y') &&
		  (buffer[1] == 'e' || buffer[1] == 'E') &&
		  (buffer[2] == 's' || buffer[2] == 'S')) 
		{
		arrayInt[i++] = 1;
		}
	      else
		{
		arrayInt[i++] = 0;
		}
	      buffer = strtok( NULL, "\t\n " );
	      }
	    } while ( i < arraySize && TYPEDSTREAM_VALUE == this->ReadLineToken() );
	  if ( i < arraySize ) 
	    {
	    Status = TYPEDSTREAM_ERROR_ARG;
	    return TYPEDSTREAM_ERROR;
	    }
	  break;
	  }
	  case TYPEDSTREAM_TYPE_BINARYBOOL: 
	  {
	  byte *arrayInt = static_cast<byte*>( array ); 
	  do 
	    {
	    buffer = strtok( BufferValue, "\t\n " );
	    int idx = 0;
	    while ( (buffer[idx]=='0') || (buffer[idx]=='1') ) 
	      {
	      if ( i >= arraySize )
		return TYPEDSTREAM_OK;
	      if ( buffer[idx] == '0' )
		arrayInt[i/8] &= ~(1<<(i%8));
	      else
		arrayInt[i/8] |= (1<<(i%8));
	      ++i;
	      ++idx;
	      }
	    } while ( i < arraySize && TYPEDSTREAM_VALUE == this->ReadLineToken() );
	  if ( i < arraySize ) 
	    {
	    Status = TYPEDSTREAM_ERROR_ARG;
	    return TYPEDSTREAM_ERROR;
	    }
	  break;
	  }
	  case TYPEDSTREAM_TYPE_FLOAT: 
	  {
	  float *arrayFloat = static_cast<float*>( array );
	  do
	    {
	    buffer = strtok( BufferValue, "\t\n " );
	    while ( buffer ) 
	      {
	      if ( *buffer == '\"' ) 
		{
		Status = TYPEDSTREAM_ERROR_TYPE;
		return TYPEDSTREAM_ERROR;
		}
	      if ( i >= arraySize )
		return TYPEDSTREAM_OK;
	      arrayFloat[i++] = static_cast<float>( atof( buffer ) );
	      buffer = strtok( NULL, "\t\n " );
	      }
	    } while ( i < arraySize && TYPEDSTREAM_VALUE == this->ReadLineToken() );
	  if ( i < arraySize ) 
	    {
	    Status = TYPEDSTREAM_ERROR_ARG;
	    return TYPEDSTREAM_ERROR;
	    }
	  break;
	  }
	case TYPEDSTREAM_TYPE_DOUBLE: 
	{
	double *arrayDouble = static_cast<double*>( array );
	do
	  {
	  buffer = strtok( BufferValue, "\t\n " );
	  while ( buffer ) 
	    {
	    if ( *buffer == '\"' ) 
	      {
	      Status = TYPEDSTREAM_ERROR_TYPE;
	      return TYPEDSTREAM_ERROR;
	      }
	    if ( i >= arraySize )
	      return TYPEDSTREAM_OK;
	    arrayDouble[i++] = atof( buffer );
	    buffer = strtok( NULL, "\t\n " );
	    }
	  } while ( i < arraySize && TYPEDSTREAM_VALUE == this->ReadLineToken() );
	if ( i < arraySize ) 
	  {
	  Status = TYPEDSTREAM_ERROR_ARG;
	  return TYPEDSTREAM_ERROR;
	  }
	break;
	}
	  case TYPEDSTREAM_TYPE_STRING: 
	  {
	  char **arrayString = static_cast<char**>( array );
	  do 
	    {
	    buffer = this->StringSplit( BufferValue );
	    while ( buffer ) 
	      {
	      if ( i >= arraySize )
		return TYPEDSTREAM_OK;
	      if ( *buffer != '\"' ) 
		{
		Status = TYPEDSTREAM_ERROR_TYPE;
		return TYPEDSTREAM_ERROR;
		}
	      char *b;
	      if (!(b = (char *)malloc(strlen(buffer)))) 
		{
		for (--i; i >= 0; i--)
		  free( arrayString[i] );
		Status = TYPEDSTREAM_ERROR_SYSTEM;
		return TYPEDSTREAM_ERROR;
		}
	      arrayString[i++] = b;
	      buffer++;
	      for (; *buffer && *buffer != '\n'; buffer++) 
		{
		if (*buffer == '\"')
		  continue;
		if (*buffer == '\\') 
		  {
		  buffer++;
		  if (*buffer == 'n') 
		    {
		    *b++ = '\n';
		    continue;
		    }
		  if (*buffer == 't') 
		    {
		    *b++ = '\t';
		    continue;
		    }
		  if (*buffer == '\0')
		    break;
		  }
		*b++ = *buffer;
		}
	      *b = '\0';
	      buffer = this->StringSplit( NULL );
	      }
	    } while ( i < arraySize && TYPEDSTREAM_VALUE == this->ReadLineToken() );
	  if ( i < arraySize ) 
	    {
	    for (--i; i >= 0; i--)
	      free( arrayString[i] );
	    Status = TYPEDSTREAM_ERROR_ARG;
	    return TYPEDSTREAM_ERROR;
	    }
	  break;
	  }
	  }
	return TYPEDSTREAM_OK;
	}
      continue;
      }
    if ( line == TYPEDSTREAM_BEGIN ) 
      {
      if ( GzFile )
	LevelStack.push( gztell( GzFile ) );
      else
	LevelStack.push( ftell( File ) );
      continue;
      }
    if ( line == TYPEDSTREAM_END ) 
      {
      if ( currentLevel == LevelStack.size() ) 
	{
	Status = TYPEDSTREAM_ERROR_NONE;
	return TYPEDSTREAM_ERROR;
	}
      LevelStack.pop();
      continue;
      }
    }
  
  return TYPEDSTREAM_ERROR;
}

TypedStreamToken
TypedStream
::ReadLineToken()
{
  if ( GzFile )
    {
    if (! gzgets( GzFile, Buffer, TYPEDSTREAM_LIMIT_BUFFER ) )
      return TYPEDSTREAM_EOF;
    }
  else
    if (! fgets( Buffer, TYPEDSTREAM_LIMIT_BUFFER, File ) )
      return TYPEDSTREAM_EOF;
  
  char* buffer;
  for ( buffer = Buffer; *buffer; buffer++)
    if (*buffer != ' ' && *buffer != '\t')
      break;

  if (*buffer == '\n' || *buffer == '!' || *buffer == '#')
    return TYPEDSTREAM_COMMENT;

  if (*buffer == '}')
    return TYPEDSTREAM_END;
  
  if (*buffer == '\"' || *buffer == '-' || *buffer == '.' || (*buffer >= '0' && *buffer <= '9')) 
    {
    BufferValue = buffer;
    return TYPEDSTREAM_VALUE;
    }
  
  if (*buffer == '_' || (*buffer >= 'a' && *buffer <= 'z') || (*buffer >= 'A' && *buffer <= 'Z')) 
    {
    BufferKey = buffer;
    for (; *buffer; buffer++)
      if (*buffer == ' ' || *buffer == '\t')
	break;
    for (; *buffer; buffer++)
      if (*buffer != ' ' && *buffer != '\t')
	break;
    BufferValue = buffer;
    if (*buffer == '{')
      return TYPEDSTREAM_BEGIN;
    
    return TYPEDSTREAM_KEY;
    }
  
  return TYPEDSTREAM_COMMENT;
}

int
TypedStream
::StringCmp( const char *s1, const char *s2 )
{
  for (; *s1 && *s2; s1++, s2++) 
    {
    if (*s1 == ' ' || *s1 == '\t' || *s1 == '\n' || *s2 == ' ' || *s2 == '\t' || *s2 == '\n') 
      {
      break;
      }
    if (*s1 == *s2)
      continue;
    if (*s1 >= 'a' && *s1 <= 'z') 
      {
      if (*s1 - ('a'-'A') == *s2)
	continue;
      }
    if (*s2 >= 'a' && *s2 <= 'z') 
      {
      if (*s2 - ('a'-'A') == *s1)
	continue;
      }
    return 1;
    }
  
  if ((*s1 == ' ' || *s1 == '\0' || *s1 == '\t' || *s1 == '\n') && (*s2 == ' ' || *s2 == '\0' || *s2 == '\t' || *s2 == '\n')) 
    {
    return 0;
    }
  
  return 1;
}

char*
TypedStream
::StringSplit( char * s1 ) const
{
  if (s1)
    SplitPosition = s1-1;
  if (SplitPosition == NULL)
    return NULL;
  
  /* skip over leading white space */
  for ( SplitPosition++; *SplitPosition == '\0' || *SplitPosition == ' ' || 
	  *SplitPosition == '\t' || *SplitPosition == '\n'; SplitPosition++ )
    if ( *SplitPosition == '\0' )
      return NULL;
  
  s1 = SplitPosition;
  
  /* find token's end */
  if ( *SplitPosition == '\"' ) 
    {
    /* skip over the special string token */
    for ( SplitPosition++; *SplitPosition && *SplitPosition != '\n' && *SplitPosition != '\t'; SplitPosition++) 
      {
      if ( *SplitPosition == '\\' && *(SplitPosition+1) ) 
	{
	SplitPosition++;
	continue;
	}
      if ( *SplitPosition == '\"' ) 
	{
	SplitPosition++;
	break;
	}
      }
    } 
  else
    {
    /* skip over a numeric value */
    for ( ; *SplitPosition; SplitPosition++ ) 
      {
      if ( *SplitPosition == ' ' || *SplitPosition == '\t' || *SplitPosition == '\n')
	break;
      }
    }
  
  if ( *SplitPosition ) 
    {
    *SplitPosition = '\0';
    } 
  else
    {
    SplitPosition = NULL;
    }
  
  return s1;
}

void
TypedStream
::DebugOutput( const char* format, ... )
{
  if ( DebugFlag != TYPEDSTREAM_DEBUG_ON ) return;

  static char buffer[1024];

  va_list args;
  va_start(args, format);
  vsnprintf( buffer, sizeof( buffer ), format, args );
  va_end(args);

  fputs( buffer, stderr );
  fputs( "\n", stderr );
}

} // namespace cmtk
