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

#include <cmtkStudy.h>
#include <cmtkFileFormat.h>
#include <cmtkVolumeIO.h>
#include <cmtkClassStream.h>

#include <string>

#ifdef HAVE_LIMITS_H
#  include <limits.h>
#endif

namespace
cmtk
{

/** \addtogroup IO */
//@{

Study::Study()
  : FileSystemPath( NULL ),
    Name( NULL ),
    Description( NULL ),
    Modality( NULL ),
    Volume( NULL ),
    CustomCalibration( false ),
    MinimumValue( 0.0 ),
    MaximumValue( 0.0 ),
    Padding( false ),
    HaveUserColorMap( false ),
    Black( 0.0 ),
    White( 0.0 )
{
}

Study::Study( const char* fileSystemPath, const char *name )
  : FileSystemPath( NULL ),
    Name( NULL ),
    Description( NULL ),
    Modality( NULL ),
    Volume( NULL ),
    CustomCalibration( false ),
    MinimumValue( 0.0 ),
    MaximumValue( 0.0 ),
    Padding( false ),
    HaveUserColorMap( false ),
    Black( 0.0 ),
    White( 0.0 )
{
  if ( fileSystemPath ) {
    FileSystemPath = strdup( fileSystemPath );
    Description = strdup( FileFormat::Describe( FileSystemPath ) );

    // cut trailing '/'s off the study path.
    char *endp = FileSystemPath + strlen( FileSystemPath ) - 1;
    while ( endp > FileSystemPath ) {
      if ( *endp == '/' ) 
	*endp = 0;
      else
	break;
    }

    this->SetMakeName( name );
  }
}

void
Study::UpdateFromVolume()
{
  const TypedArray *dataArray = this->Volume->GetData();
  if ( dataArray ) 
    {
    if ( dataArray->GetRange( MinimumValue, MaximumValue ) )
      {
      const Types::DataItem perc01 = dataArray->GetPercentile( 0.01, 1024 );
      const Types::DataItem perc99 = dataArray->GetPercentile( 0.99, 1024 );

      Black = std::min( std::max( Black, perc01 ), MaximumValue );
      White = std::max( std::min( White, perc99 ), MinimumValue );
      }  
    else
      {
      Black = White = MinimumValue = MaximumValue = 0.0;
      }
    }
}
  
const char*
Study::SetMakeName( const char* name, const int suffix )
{
  if ( name )
    {
    if ( suffix )
      {
      char fullname[PATH_MAX];
      snprintf( fullname, sizeof( fullname ), "%s <%d>", name, suffix );
      this->SetName( fullname );
      }
    else
      {
      this->SetName( name );
      }
    }
  else
    {
    char buffer[PATH_MAX];
    strncpy( buffer, FileSystemPath, PATH_MAX );
    
    char* lastChar = buffer + strlen( buffer ) - 1;
    while ( (lastChar != buffer) && // not yet empty string
	    (*lastChar =='/') ) 
      {
      *lastChar = 0;
      --lastChar;
      }
    
    const char *slash = strrchr( buffer, '/' );
    if ( slash )
      strcpy( buffer, slash+1 );
    else 
      strcpy( buffer, FileSystemPath );
    
    char* dot = strchr( buffer, '.' );
    if ( dot )
      *dot = 0;
    else
      dot = buffer + strlen(buffer);
    
    if ( suffix ) 
      {
      snprintf( dot, sizeof( buffer ) - (dot-buffer), "<%d>", suffix );
      }

    this->SetName( buffer );
    }

  return Name;
}

Study::~Study() 
{
  if ( FileSystemPath ) free( FileSystemPath );
  if ( Description ) free( Description );
  if ( Name ) free( Name );
}

bool 
Study::ReadVolume( const bool reRead, const char* orientation )
{
  UniformVolume::SmartPtr oldVolume( NULL );

  if ( Volume && reRead ) 
    {
    oldVolume = Volume;
    Volume = UniformVolume::SmartPtr( NULL );
    }
  
  if ( !Volume ) 
    {
    if ( orientation )
      Volume = UniformVolume::SmartPtr( VolumeIO::ReadOriented( FileSystemPath, orientation ) );
    else
      Volume = UniformVolume::SmartPtr( VolumeIO::Read( FileSystemPath ) );
    
    if ( Volume ) 
      {
      this->SetDims( Volume->GetDims( AXIS_X ), Volume->GetDims( AXIS_Y ), Volume->GetDims( AXIS_Z ) );
      DisplayedImageIndex = Volume->GetDims( AXIS_Z ) / 2 ;
      ZoomFactor = 1;
      const TypedArray *dataArray = Volume->GetData();
      if ( dataArray ) 
	{
	if ( dataArray->GetRange( MinimumValue, MaximumValue ) )
	  {
	  Black = dataArray->GetPercentile( 0.01, 1024 );
	  White = dataArray->GetPercentile( 0.99, 1024 );
	  }
	else
	  {
	  Black = White = MinimumValue = MaximumValue = 0.0;
	  }
	StandardColormap = 0;
	ReverseColormap = false;
	}
      }
    }
  
  if ( Volume ) 
    {
    if ( this->m_LandmarkList ) 
      {
      Volume->m_LandmarkList = this->m_LandmarkList;
      }
    if ( Volume->GetData() ) 
      {
      return true;
      }
    }
  
  Volume = oldVolume;
  return false;
}

Study* 
Study::Read( const char* path )
{
  Study* study = NULL;
  if ( FileFormat::Identify( path ) == FILEFORMAT_STUDY ) 
    {    
    ClassStream stream( path, "images", ClassStream::READ );
    if ( stream.IsValid() ) 
      {
      stream >> study;
      }
    stream.Close();
    
    if ( study ) 
      {
      study->SetFileSystemPath( path );
      study->SetMakeName(); // create a default name for this study
      }
    
    stream.Open( path, "landmarks", ClassStream::READ );
    if ( stream.IsValid() )
      {
      LandmarkList::SmartPtr ll( NULL );
      stream >> ll;
      if ( ll ) 
	{
	study->m_LandmarkList = ll;
	}
      }
    stream.Close();
    } 
  else
    {
    study = new Study( path );
    }
  
  return study;
}

bool
Study::Write() const
{
  ClassStream stream( FileSystemPath, "images", ClassStream::WRITE );
  if ( ! stream.IsValid() ) return false;

  stream << *this;
  stream.Close();

  stream.Open( FileSystemPath, "landmarks", ClassStream::WRITE );
  if ( ! stream.IsValid() ) return false;

  stream << this->m_LandmarkList;
  stream.Close();

  return true;
}

bool
Study::WriteTo( const char* path )
{
  ClassStream stream( path, "images", ClassStream::WRITE );
  if ( ! stream.IsValid() ) return false;

  stream << *this;
  this->SetFileSystemPath( path );
  stream.Close();

  if ( this->m_LandmarkList ) 
    {
    stream.Open( FileSystemPath, "landmarks", ClassStream::WRITE );
    if ( ! stream.IsValid() ) return false;
    
    stream << this->m_LandmarkList;
    stream.Close();
    }

  return true;
}

void
Study::CopyColormap( const Study* other )
{
  MinimumValue = other->MinimumValue;
  MaximumValue = other->MaximumValue;
  StandardColormap = other->StandardColormap;
  ReverseColormap = other->ReverseColormap;
  Black = other->Black;
  White = other->White;
  Gamma = other->Gamma;
}

} // namespace cmtk
