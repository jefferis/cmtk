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

#include "cmtkStudy.h"

#include <IO/cmtkFileFormat.h>
#include <IO/cmtkVolumeIO.h>

#include <string>
#include <limits.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

Study::Study()
  : m_Modality( NULL ),
    m_Volume( NULL ),
    m_MinimumValue( 0.0 ),
    m_MaximumValue( 0.0 ),
    m_Padding( false ),
    m_PaddingValue( 0.0 ),
    m_HaveUserColorMap( false ),
    m_StandardColormap( 0 ),
    m_ReverseColormap( false ),
    m_Black( 0.0 ),
    m_White( 0.0 ),
    m_Gamma( 1.0 ),
    m_DisplayedImageIndex( -1 ),
    m_ZoomFactor( 1 ),
    m_SliceNormal( 2 )
{
}

Study::Study( const std::string& fileSystemPath, const std::string& name )
  : m_Modality( NULL ),
    m_Volume( NULL ),
    m_MinimumValue( 0.0 ),
    m_MaximumValue( 0.0 ),
    m_Padding( false ),
    m_PaddingValue( 0.0 ),
    m_HaveUserColorMap( false ),
    m_StandardColormap( 0 ),
    m_ReverseColormap( false ),
    m_Black( 0.0 ),
    m_White( 0.0 ),
    m_Gamma( 1.0 ),
    m_DisplayedImageIndex( -1 ),
    m_ZoomFactor( 1 ),
    m_SliceNormal( 2 )
{
  if ( ! fileSystemPath.empty() ) 
    {
    this->m_FileSystemPath = fileSystemPath;
    this->m_Description = FileFormat::Describe( this->m_FileSystemPath );

    // cut trailing '/'s off the study path.
    const size_t lastChar = this->m_FileSystemPath.find_last_not_of( "/" );
    if ( lastChar != std::string::npos )
      {
      this->m_FileSystemPath = this->m_FileSystemPath.substr( 0, lastChar+1 );
      }
    
    this->SetMakeName( name );
  }
}

std::string
Study::SetMakeName( const std::string& name, const int suffix )
{
  char suffixStr[10];
  snprintf( suffixStr, 9, "<%d>", suffix );
  if ( !name.empty() )
    {
    if ( suffix )
      {
      this->SetName( name + suffixStr );
      }
    else
      {
      this->SetName( name );
      }
    }
  else
    {
    std::string studyName = name;

    const size_t lastChar = studyName.find_last_not_of( "/" );
    if ( lastChar != std::string::npos )
      {
      studyName = studyName.substr( 0, lastChar+1 );
      }

    const size_t lastSlash = studyName.rfind( "/" );
    if ( lastSlash != std::string::npos )
      {
      studyName = studyName.substr( lastSlash+1 );
      }
    else
      {
      studyName = this->m_FileSystemPath;
      }
    
    const size_t dot = studyName.find( "." );
    if ( dot != std::string::npos )
      {
      studyName = studyName.substr( 0, dot );
      }

    if ( suffix )
      studyName = studyName + suffixStr;

    this->SetName( studyName );
    }

  return this->m_Name;
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
	const Types::DataItemRange range = dataArray->GetRange();
	this->m_MinimumValue = range.m_LowerBound;
	this->m_MaximumValue = range.m_UpperBound;

	this->m_Black = dataArray->GetPercentile( 0.01, 1024 );
	this->m_White = dataArray->GetPercentile( 0.99, 1024 );

	this->m_StandardColormap = 0;
	this->m_ReverseColormap = false;
	}
      }
    }
  
  if ( this->m_Volume && this->m_Volume->GetData() ) 
    {
    return true;
    }
  
  this->m_Volume = oldVolume;
  return false;
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
