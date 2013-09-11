/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012, 2013 SRI International
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

#include "cmtkTypedStreamInput.h"

#include <System/cmtkFileUtils.h>
#include <System/cmtkConsole.h>

#include <string.h>
#include <stdlib.h>
#include <limits.h>

#ifdef HAVE_MALLOC_H
#  include <malloc.h>
#endif

namespace
cmtk
{

TypedStreamInput::TypedStreamInput
( const std::string& filename )
{
  this->Open( filename );
}

TypedStreamInput::TypedStreamInput
( const std::string& dir, const std::string& archive )
{
  this->Open( dir, archive );
}

TypedStreamInput
::~TypedStreamInput()
{
  this->Close();
}

void 
TypedStreamInput
::Open
( const std::string& dir, const std::string& archive )
{
  static char fname[PATH_MAX];
  
  // If "dir" parameter is empty, use current directory instead.
  if ( dir != "" ) 
    {
    if ( static_cast<size_t>( snprintf( fname, sizeof( fname ), "%s%c%s", dir.c_str(), CMTK_PATH_SEPARATOR, archive.c_str() ) ) >= sizeof( fname ) )
      {
      StdErr << "WARNING: length of path exceeds system PATH_MAX in TypedStreamInput::Open and will be truncated.\n";
      }
    } 
  else 
    {
    if ( static_cast<size_t>( snprintf( fname, sizeof( fname ), "%s", archive.c_str() ) ) >= sizeof( fname ) )
      {
      StdErr << "WARNING: length of path exceeds system PATH_MAX in TypedStreamInput::Open and will be truncated.\n";
      }
    }
  
  this->Open( fname );
}

void 
TypedStreamInput
::Open
( const std::string& filename )
{
  this->m_Status = Self::ERROR_NONE;
  this->Close();
  
#ifdef _MSC_VER
  const char *modestr = "rb";
#else
  const char *modestr = "r";
#endif

  if ( ! ( File = fopen( filename.c_str(), modestr ) ) ) 
    {
    const std::string gzName = filename + ".gz";
    GzFile = gzopen( gzName.c_str(), modestr );
    if ( ! GzFile ) 
      {
      StdErr << "ERROR: could not open file \"" << filename << "\" with mode \"" << modestr << "\"\n";
      this->m_Status = Self::ERROR_SYSTEM;
      return;
      }
    }
  
  if ( GzFile ) 
    {
    if (! gzgets( GzFile, Buffer, sizeof( this->Buffer ) ) ) 
      {
      this->m_Status = Self::ERROR_FORMAT;
      gzclose( GzFile );
      return;
      }
    } 
  else
    if (! fgets( Buffer, sizeof( this->Buffer) , File ) ) 
      {
      this->m_Status = Self::ERROR_FORMAT;
      fclose( File );
      File = NULL;
      return;
      }
  
  if (Buffer[0] != '!' && Buffer[0] != '#') 
    {
    this->m_Status = Self::ERROR_FORMAT;
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
    return;
    }
  int releaseMinor, releaseMajor;
  if (2 != sscanf( Buffer+1, " TYPEDSTREAM %d.%d", &releaseMajor, &releaseMinor)) 
    {
    this->m_Status = Self::ERROR_FORMAT;
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
}

void
TypedStreamInput
::Close()
{
  if ( File || GzFile )
    {
    while ( ! LevelStack.empty() ) 
      {
      LevelStack.pop();      
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
  
  this->m_Status = Self::ERROR_NONE;
  SplitPosition = NULL;
}

TypedStreamInput::Condition 
TypedStreamInput
::Begin()
{
  if ( !File && !GzFile)
    {
    this->m_Status = Self::ERROR_INVALID;
    return Self::CONDITION_ERROR;
    }

  if ( GzFile )
    gzseek( GzFile, LevelStack.top(), SEEK_SET );
  else
    fseek( File, LevelStack.top(), SEEK_SET );
  return Self::CONDITION_OK;
}

TypedStreamInput::Condition
TypedStreamInput
::End()
{
  if ( ! File && ! GzFile )
    {
    this->m_Status = Self::ERROR_INVALID;
    return Self::CONDITION_ERROR;
    }
  
  if ( LevelStack.empty() ) 
    {
    // end without begin
    this->m_Status = Self::ERROR_LEVEL;
    return Self::CONDITION_ERROR;
    }
  
  Self::Token token;
  int currentLevel = 1;
  while ( currentLevel && (Self::TOKEN_EOF != ( token = this->ReadLineToken() ) ) ) 
    {
    if ( token == Self::TOKEN_BEGIN ) 
      {
      this->DebugOutput( "Skipping section %s at level %d.", BufferKey, currentLevel );
      ++currentLevel;
      }
    else
      if ( token == Self::TOKEN_END ) 
	{
	this->DebugOutput( "Leaving section %d.", currentLevel );
	--currentLevel;
	}
      }
  
  LevelStack.pop();
  
  return Self::CONDITION_OK;
}
  
TypedStreamInput::Condition
TypedStreamInput
::Seek
( const char* section, const bool forward )
{
  if ( ! File && ! GzFile )
    {
    this->m_Status = Self::ERROR_INVALID;
    return Self::CONDITION_ERROR;
    }

  if ( ! section ) 
    {
    this->m_Status = Self::ERROR_ARG;
    return Self::CONDITION_ERROR;
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
  Self::Token token;
  while ( Self::TOKEN_EOF != ( token = this->ReadLineToken() ) ) 
    {
    if ( token == Self::TOKEN_BEGIN ) 
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
	  return Self::CONDITION_OK;
	  }
	if ( currentLevel == LevelStack.size()-1 ) 
	  {
	  LevelStack.pop();
	  if ( GzFile )
	    LevelStack.push( gztell( GzFile ) );
	  else
	    LevelStack.push( ftell( File ) );
	  return Self::CONDITION_OK;
	  }
	}
      ++ currentLevel;
      }
    if ( token == Self::TOKEN_END ) 
      {
      this->DebugOutput( "Leaving section %d.", currentLevel );
      if ( ! currentLevel ) 
	{
	this->m_Status = Self::ERROR_LEVEL;
	return Self::CONDITION_ERROR;
	}
      if ( currentLevel < initialLevel ) 
	{
	this->m_Status = Self::ERROR_NONE;
	return Self::CONDITION_ERROR;
	}
      -- currentLevel;
      }
    }
  
  this->DebugOutput( "Section %s not found.", section );
  this->m_Status = Self::ERROR_NONE;
  return Self::CONDITION_ERROR;
}

TypedStreamInput::Condition
TypedStreamInput
::Rewind()
{
  if ( ! File && ! GzFile )
    {
    this->m_Status = Self::ERROR_INVALID;
    return Self::CONDITION_ERROR;
    }
  
  if ( !LevelStack.empty() )
    LevelStack.pop();
  
  if ( LevelStack.empty() ) 
    {
    if ( GzFile ) 
      {
      if ( -1 == gzseek( GzFile, 0, SEEK_SET ) ) 
	{
	this->m_Status = Self::ERROR_SYSTEM;
	return Self::CONDITION_ERROR;
	}
      } 
    else
      if ( -1 == fseek( File, 0, SEEK_SET ) )
	{
	this->m_Status = Self::ERROR_SYSTEM;
	return Self::CONDITION_ERROR;
	}
    } 
  else 
    {
    if ( GzFile ) 
      {
      if (-1 == gzseek( GzFile, LevelStack.top(), SEEK_SET)) 
	{
	this->m_Status = Self::ERROR_SYSTEM;
	return Self::CONDITION_ERROR;
	}
      }
    else
      if (-1 == fseek( File, LevelStack.top(), SEEK_SET)) 
	{
	this->m_Status = Self::ERROR_SYSTEM;
	return Self::CONDITION_ERROR;
	}
    }
  return Self::CONDITION_OK;
}

bool
TypedStreamInput
::ReadBool( const char* key, const bool defaultValue, const bool forward )
{
  int value;
  
  if ( this->GenericReadArray( key, Self::TYPE_BOOL, &value, 1, forward ) != Self::CONDITION_OK )
    if ( this->GenericReadArray( key, Self::TYPE_INT, &value, 1, forward ) != Self::CONDITION_OK )
      return defaultValue;
  
  return (value != 0);
}

TypedStreamInput::Condition 
TypedStreamInput
::ReadBoolArray( const char* key, byte *const array, const int size, const bool forward )
{
  return this->GenericReadArray( key, Self::TYPE_BINARYBOOL, array, size, forward );
}

int 
TypedStreamInput
::ReadInt( const char* key, const int defaultValue, const bool forward )
{
  int value = defaultValue;
  if ( this->GenericReadArray( key, Self::TYPE_INT, &value, 1, forward ) != Self::CONDITION_OK )
    return defaultValue;
  
  return value;
}

TypedStreamInput::Condition 
TypedStreamInput
::ReadIntArray( const char* key, int *const array, const int size, const bool forward )
{
  return this->GenericReadArray( key, Self::TYPE_INT, array, size, forward );
}

float
TypedStreamInput
::ReadFloat( const char* key, const float defaultValue, const bool forward )
{
  float value = defaultValue;
  if ( this->GenericReadArray( key, Self::TYPE_FLOAT, &value, 1, forward ) != Self::CONDITION_OK )
    return defaultValue;
  
  return value;
}

TypedStreamInput::Condition
TypedStreamInput
::ReadFloatArray( const char* key, float *const array, const int size, const bool forward )
{
  return this->GenericReadArray( key, Self::TYPE_FLOAT, array, size, forward );
}

double
TypedStreamInput
::ReadDouble( const char* key, const double defaultValue, const bool forward )
{
  double value = defaultValue;
  if ( this->GenericReadArray( key, Self::TYPE_DOUBLE, &value, 1, forward ) != Self::CONDITION_OK )
    return defaultValue;
  
  return value;
}

TypedStreamInput::Condition
TypedStreamInput
::ReadDoubleArray( const char* key, double *const array, const int size, const bool forward )
{
  return this->GenericReadArray( key, Self::TYPE_DOUBLE, array, size, forward );
}

char*
TypedStreamInput
::ReadString( const char* key, const char *defaultValue, const bool forward )
{
  char *value;
  if ( this->GenericReadArray( key, Self::TYPE_STRING, &value, 1, forward ) != Self::CONDITION_OK )
    {
    if ( defaultValue )
      return strdup( defaultValue );
    else
      return NULL;
    }
  
  return value;
}

TypedStreamInput::Condition
TypedStreamInput
::GenericReadArray( const char * key, const int type, void *const array, const int arraySize, const bool forward )
{
  if (!array || arraySize < 1) 
    {
    this->m_Status = Self::ERROR_ARG;
    return Self::CONDITION_ERROR;
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
  while ( Self::TOKEN_EOF != ( line = this->ReadLineToken() ) ) 
    {
    if ( line == Self::TOKEN_KEY) 
      {
      if ( (currentLevel == LevelStack.size()) && (this->StringCmp( BufferKey, key )) == 0 ) 
	{
	int i = 0;
	switch (type) 
	  {
	  case Self::TYPE_INT: 
	  {
	  int *arrayInt = static_cast<int*>( array ); 
	  do
	    {
	    buffer = strtok( BufferValue, "\t\n " );
	    while ( buffer ) 
	      {
	      if ( *buffer == '\"' ) 
		{
		this->m_Status = Self::ERROR_TYPE;
		return Self::CONDITION_ERROR;
		}
	      if ( i >= arraySize )
		return Self::CONDITION_OK;
	      arrayInt[i++] = atoi(buffer);
	      buffer = strtok( NULL, "\t\n " );
	      }
	    } while ( i < arraySize && Self::TOKEN_VALUE == this->ReadLineToken() );
	  if ( i < arraySize) 
	    {
	    this->m_Status = Self::ERROR_ARG;
	    return Self::CONDITION_ERROR;
	    }
	  break;
	  }
	  case Self::TYPE_BOOL: 
	  {
	  int *arrayInt = static_cast<int*>( array ); 
	  if (arraySize != 1) 
	    {
	    this->m_Status = Self::ERROR_ARG;
	    return Self::CONDITION_ERROR;
	    }
	  do 
	    {
	    buffer = strtok( BufferValue, "\t\n " );
	    while ( buffer ) 
	      {
	      if ( i >= arraySize )
		return Self::CONDITION_OK;
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
	    } while ( i < arraySize && Self::TOKEN_VALUE == this->ReadLineToken() );
	  if ( i < arraySize ) 
	    {
	    this->m_Status = Self::ERROR_ARG;
	    return Self::CONDITION_ERROR;
	    }
	  break;
	  }
	  case Self::TYPE_BINARYBOOL: 
	  {
	  byte *arrayInt = static_cast<byte*>( array ); 
	  do 
	    {
	    buffer = strtok( BufferValue, "\t\n " );
	    int idx = 0;
	    while ( (buffer[idx]=='0') || (buffer[idx]=='1') ) 
	      {
	      if ( i >= arraySize )
		return Self::CONDITION_OK;
	      if ( buffer[idx] == '0' )
		arrayInt[i/8] &= ~(1<<(i%8));
	      else
		arrayInt[i/8] |= (1<<(i%8));
	      ++i;
	      ++idx;
	      }
	    } while ( i < arraySize && Self::TOKEN_VALUE == this->ReadLineToken() );
	  if ( i < arraySize ) 
	    {
	    this->m_Status = Self::ERROR_ARG;
	    return Self::CONDITION_ERROR;
	    }
	  break;
	  }
	  case Self::TYPE_FLOAT: 
	  {
	  float *arrayFloat = static_cast<float*>( array );
	  do
	    {
	    buffer = strtok( BufferValue, "\t\n " );
	    while ( buffer ) 
	      {
	      if ( *buffer == '\"' ) 
		{
		this->m_Status = Self::ERROR_TYPE;
		return Self::CONDITION_ERROR;
		}
	      if ( i >= arraySize )
		return Self::CONDITION_OK;
	      arrayFloat[i++] = static_cast<float>( atof( buffer ) );
	      buffer = strtok( NULL, "\t\n " );
	      }
	    } while ( i < arraySize && Self::TOKEN_VALUE == this->ReadLineToken() );
	  if ( i < arraySize ) 
	    {
	    this->m_Status = Self::ERROR_ARG;
	    return Self::CONDITION_ERROR;
	    }
	  break;
	  }
	case Self::TYPE_DOUBLE: 
	{
	double *arrayDouble = static_cast<double*>( array );
	do
	  {
	  buffer = strtok( BufferValue, "\t\n " );
	  while ( buffer ) 
	    {
	    if ( *buffer == '\"' ) 
	      {
	      this->m_Status = Self::ERROR_TYPE;
	      return Self::CONDITION_ERROR;
	      }
	    if ( i >= arraySize )
	      return Self::CONDITION_OK;
	    arrayDouble[i++] = atof( buffer );
	    buffer = strtok( NULL, "\t\n " );
	    }
	  } while ( i < arraySize && Self::TOKEN_VALUE == this->ReadLineToken() );
	if ( i < arraySize ) 
	  {
	  this->m_Status = Self::ERROR_ARG;
	  return Self::CONDITION_ERROR;
	  }
	break;
	}
	  case Self::TYPE_STRING: 
	  {
	  char **arrayString = static_cast<char**>( array );
	  do 
	    {
	    buffer = this->StringSplit( BufferValue );
	    while ( buffer ) 
	      {
	      if ( i >= arraySize )
		return Self::CONDITION_OK;
	      if ( *buffer != '\"' ) 
		{
		this->m_Status = Self::ERROR_TYPE;
		return Self::CONDITION_ERROR;
		}
	      char *b;
	      if (!(b = (char *)malloc(1+strlen(buffer)))) 
		{
		for (--i; i >= 0; i--)
		  free( arrayString[i] );
		this->m_Status = Self::ERROR_SYSTEM;
		return Self::CONDITION_ERROR;
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
	    } while ( i < arraySize && Self::TOKEN_VALUE == this->ReadLineToken() );
	  if ( i < arraySize ) 
	    {
	    for (--i; i >= 0; i--)
	      free( arrayString[i] );
	    this->m_Status = Self::ERROR_ARG;
	    return Self::CONDITION_ERROR;
	    }
	  break;
	  }
	  }
	return Self::CONDITION_OK;
	}
      continue;
      }
    if ( line == Self::TOKEN_BEGIN ) 
      {
      if ( GzFile )
	LevelStack.push( gztell( GzFile ) );
      else
	LevelStack.push( ftell( File ) );
      continue;
      }
    if ( line == Self::TOKEN_END ) 
      {
      if ( currentLevel == LevelStack.size() ) 
	{
	this->m_Status = Self::ERROR_NONE;
	return Self::CONDITION_ERROR;
	}
      LevelStack.pop();
      continue;
      }
    }
  
  return Self::CONDITION_ERROR;
}

TypedStreamInput::Token
TypedStreamInput
::ReadLineToken()
{
  if ( GzFile )
    {
    if (! gzgets( GzFile, Buffer, sizeof( this->Buffer ) ) )
      return Self::TOKEN_EOF;
    }
  else
    if (! fgets( Buffer, sizeof( this->Buffer ), File ) )
      return Self::TOKEN_EOF;
  
  char* buffer;
  for ( buffer = Buffer; *buffer; buffer++)
    if (*buffer != ' ' && *buffer != '\t')
      break;

  if (*buffer == '\n' || *buffer == '!' || *buffer == '#')
    return Self::TOKEN_COMMENT;

  if (*buffer == '}')
    return Self::TOKEN_END;
  
  if (*buffer == '\"' || *buffer == '-' || *buffer == '.' || (*buffer >= '0' && *buffer <= '9')) 
    {
    BufferValue = buffer;
    return Self::TOKEN_VALUE;
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
      return Self::TOKEN_BEGIN;
    
    return Self::TOKEN_KEY;
    }
  
  return Self::TOKEN_COMMENT;
}

} // namespace cmtk
