/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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
  : m_FileSystemPath( NULL ),
    m_Name( NULL ),
    m_Description( NULL ),
    m_Modality( NULL ),
    m_Volume( NULL ),
    m_CustomCalibration( false ),
    m_MinimumValue( 0.0 ),
    m_MaximumValue( 0.0 ),
    m_Padding( false ),
    m_HaveUserColorMap( false ),
    m_Black( 0.0 ),
    m_White( 0.0 )
{
}

Study::Study( const char* fileSystemPath, const char *name )
  : m_FileSystemPath( NULL ),
    m_Name( NULL ),
    m_Description( NULL ),
    m_Modality( NULL ),
    m_Volume( NULL ),
    m_CustomCalibration( false ),
    m_MinimumValue( 0.0 ),
    m_MaximumValue( 0.0 ),
    m_Padding( false ),
    m_HaveUserColorMap( false ),
    m_Black( 0.0 ),
    m_White( 0.0 )
{
  if ( fileSystemPath ) 
    {
    this->m_FileSystemPath = strdup( fileSystemPath );
    this->m_Description = strdup( FileFormat::Describe( this->m_FileSystemPath ) );

    // cut trailing '/'s off the study path.
    char *endp = this->m_FileSystemPath + strlen( this->m_FileSystemPath ) - 1;
    while ( endp > this->m_FileSystemPath ) 
      {
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
  const TypedArray *dataArray = this->m_Volume->GetData();
  if ( dataArray ) 
    {
    if ( dataArray->GetRange( this->m_MinimumValue, this->m_MaximumValue ) )
      {
      const Types::DataItem perc01 = dataArray->GetPercentile( 0.01, 1024 );
      const Types::DataItem perc99 = dataArray->GetPercentile( 0.99, 1024 );

      this->m_Black = std::min( std::max( this->m_Black, perc01 ), this->m_MaximumValue );
      this->m_White = std::max( std::min( this->m_White, perc99 ), this->m_MinimumValue );
      }  
    else
      {
      this->m_Black = this->m_White = this->m_MinimumValue = this->m_MaximumValue = 0.0;
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
    strncpy( buffer, this->m_FileSystemPath, PATH_MAX );
    
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
      strcpy( buffer, this->m_FileSystemPath );
    
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

  return this->m_Name;
}

Study::~Study() 
{
  free( this->m_FileSystemPath );
  free( this->m_Description );
  free( this->m_Name );
}

bool 
Study::ReadVolume( const bool reRead, const char* orientation )
{
  UniformVolume::SmartPtr oldVolume( NULL );

  if ( this->m_Volume && reRead ) 
    {
    oldVolume = this->m_Volume;
    this->m_Volume = UniformVolume::SmartPtr( NULL );
    }
  
  if ( !this->m_Volume ) 
    {
    if ( orientation )
      this->m_Volume = UniformVolume::SmartPtr( VolumeIO::ReadOriented( this->m_FileSystemPath, orientation ) );
    else
      this->m_Volume = UniformVolume::SmartPtr( VolumeIO::Read( this->m_FileSystemPath ) );
    
    if ( this->m_Volume ) 
      {
      this->m_Dims = this->m_Volume->GetDims();
      this->m_DisplayedImageIndex = this->m_Dims[AXIS_Z] / 2 ;
      this->m_ZoomFactor = 1;
      const TypedArray *dataArray = this->m_Volume->GetData();
      if ( dataArray ) 
	{
	if ( dataArray->GetRange( this->m_MinimumValue, this->m_MaximumValue ) )
	  {
	  this->m_Black = dataArray->GetPercentile( 0.01, 1024 );
	  this->m_White = dataArray->GetPercentile( 0.99, 1024 );
	  }
	else
	  {
	  this->m_Black = this->m_White = this->m_MinimumValue = this->m_MaximumValue = 0.0;
	  }
	this->m_StandardColormap = 0;
	this->m_ReverseColormap = false;
	}
      }
    }
  
  if ( this->m_Volume ) 
    {
    if ( this->m_LandmarkList ) 
      {
      this->m_Volume->m_LandmarkList = this->m_LandmarkList;
      }
    if ( this->m_Volume->GetData() ) 
      {
      return true;
      }
    }
  
  this->m_Volume = oldVolume;
  return false;
}

Study* 
Study::Read( const char* path )
{
  return new Study( path );
}

void
Study::CopyColormap( const Study* other )
{
  this->m_MinimumValue = other->m_MinimumValue;
  this->m_MaximumValue = other->m_MaximumValue;
  this->m_StandardColormap = other->m_StandardColormap;
  this->m_ReverseColormap = other->m_ReverseColormap;
  this->m_Black = other->m_Black;
  this->m_White = other->m_White;
  this->m_Gamma = other->m_Gamma;
}

} // namespace cmtk
