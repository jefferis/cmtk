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

#include <cmtkClassStream.h>

#include <cmtkStudyImageSet.h>
#include <cmtkStudyImageSetRaw.h>


namespace
cmtk
{

/** \addtogroup IO */
//@{

ClassStream& 
ClassStream::operator << ( const Study *study )
{
  if ( !this->IsValid() ) return *this;

  this->Begin( "imageserie" );
    
  FileFormatID imageFormat = study->GetImageFormat();
  this->WriteString( "format", FileFormatName[ imageFormat ] );

  const StudyImageSetRaw* studyRaw = dynamic_cast<const StudyImageSetRaw*>( study );
  if ( studyRaw ) 
    {
    this->WriteInt( "offset", studyRaw->GetHeaderLength() );
    this->WriteInt( "bytesperpixel", studyRaw->GetBytesPerPixel() );
    this->WriteBool( "big_endian", studyRaw->GetBigEndian() );
    this->WriteBool( "little_endian", !studyRaw->GetBigEndian() );
    this->WriteBool( "swapbytes", !studyRaw->GetBigEndian() );
    this->WriteBool( "signed", studyRaw->GetSigned() );
    }
  
  this->WriteString( "modality", study->GetModality() );

  this->WriteInt( "width", study->GetDims( AXIS_X ) );
  this->WriteInt( "height", study->GetDims( AXIS_Y ) );
  this->WriteInt( "depth", study->GetDims( AXIS_Z ) );

  this->WriteInt( "minimum", static_cast<int>( study->GetMinimumValue() ) );
  this->WriteInt( "maximum", static_cast<int>( study->GetMaximumValue() ) );

  this->WriteBool( "padding", study->GetPadding() );
  this->WriteInt( "padding_value", static_cast<int>( study->GetPaddingValue() ) );
  
  this->WriteInt( "black", static_cast<int>( study->GetBlack() ) );
  this->WriteInt( "white", static_cast<int>( study->GetWhite() ) );
  this->WriteDouble( "gamma", study->GetGamma() );
  this->WriteBool( "reverse", study->GetReverseColormap() );
    
  this->WriteBool( "custom", study->GetCustomCalibration() );
  
  const UniformVolume* volume = study->GetVolume();
  if ( study->GetCustomCalibration() || !volume ) 
    {
    this->WriteDouble( "calibrationx", study->GetCalibration( AXIS_X ) );
    this->WriteDouble( "calibrationy", study->GetCalibration( AXIS_Y ) );
    this->WriteDouble( "slicedistance", study->GetCalibration( AXIS_Z ) );
    } 
  else 
    {
    this->WriteDouble( "calibrationx", volume->GetDelta( AXIS_X, 0 ) );
    this->WriteDouble( "calibrationy", volume->GetDelta( AXIS_Y, 0 ) );
    this->WriteDouble( "slicedistance", volume->GetDelta( AXIS_Z, 0 ) );
    }
  
  const StudyImageSet* studySet = dynamic_cast<const StudyImageSet*>( study );
  
  if ( studySet )
    this->WriteString( "imagepath", studySet->GetImageDirectory() );

  this->End();
  
  if ( studySet ) 
    {
    StudyImageSet::const_iterator it = studySet->begin();
    
    unsigned int idx = 0;
    while ( it != studySet->end() ) 
      {
      this->Begin( "image" );
      this->WriteString( "name", it->c_str() );
      this->WriteDouble( "calibrationx", volume->GetDelta( AXIS_X, 0 ) );
      this->WriteDouble( "calibrationy", volume->GetDelta( AXIS_Y, 0 ) );
      this->WriteDouble( "tablepos", volume->GetPlaneCoord( AXIS_Z, idx ) );
      this->End();
      ++it;
      ++idx;
      }
    }
  
  return *this;
}

ClassStream& 
ClassStream::operator >> ( Study*& study )
{
  if ( !this->IsValid() ) return *this;

  StudyImageSet* newStudy = NULL;
  
  if ( this->Seek( "imageserie" ) != TYPEDSTREAM_OK ) return *this;
    
  char *tmpStr = this->ReadString( "format" );
  if ( tmpStr ) 
    {
    FileFormatID imageFormat = FileFormat::GetID( tmpStr );
    switch ( imageFormat ) 
      {
      case FILEFORMAT_RAW:
      case FILEFORMAT_RAW3D: 
      {
      StudyImageSetRaw* newStudyRaw = new StudyImageSetRaw;
      newStudyRaw->SetMultiFile( imageFormat == FILEFORMAT_RAW );
      newStudyRaw->SetHeaderLength( this->ReadInt( "offset" ) );
      newStudyRaw->SetBytesPerPixel( this->ReadInt( "bytesperpixel" ) );
      newStudyRaw->SetBigEndian( this->ReadBool( "big_endian", !this->ReadBool( "swapbytes" ) ) );
      newStudyRaw->SetSigned( this->ReadBool( "signed" ) );
      newStudy = newStudyRaw;
      }
      break;
      default:
	newStudy = new StudyImageSet;
	newStudy->SetMultiFile( true );
	break;
      }
    newStudy->SetImageFormat( imageFormat );
    free( tmpStr );
    } 
  else
    {
    return *this;
    }
  
  tmpStr = this->ReadString( "imagepath" );
  if ( tmpStr )
    {
    newStudy->SetImageDirectory( tmpStr );
    free( tmpStr );
    }
  
  tmpStr = this->ReadString( "modality" );
  if ( tmpStr ) 
    {
    newStudy->SetModality( tmpStr );
    free( tmpStr );
    } 
  else
    {
    newStudy->SetModality( "UNKNOWN" );
    }
  
  unsigned int Width = this->ReadInt( "width" );
  unsigned int Height = this->ReadInt( "height" );
  unsigned int Depth = this->ReadInt( "depth", 0 );
  newStudy->SetDims( Width, Height, Depth );
    
  newStudy->SetMinimumValue( this->ReadInt( "minimum" ) );
  newStudy->SetMaximumValue( this->ReadInt( "maximum" ) );

  newStudy->SetPadding( this->ReadBool( "padding", false ) );
  newStudy->SetPaddingValue( this->ReadInt( "padding_value" ) );
    
  tmpStr = this->ReadString( "palette" );
  if ( tmpStr )
    {
    newStudy->SetStandardColormap( 0 );
    free( tmpStr );
    }
  
  newStudy->SetBlack( this->ReadInt( "black" ) );
  newStudy->SetWhite( this->ReadInt( "white" ) );
  newStudy->SetGamma( this->ReadDouble( "gamma" ) );
  newStudy->SetReverseColormap( this->ReadBool( "reverse" ) );
  
  newStudy->SetCustomCalibration( this->ReadBool( "custom" ) );
    
  newStudy->SetCalibration( this->ReadDouble( "calibrationx", 0 ), this->ReadDouble( "calibrationy", 0 ), this->ReadDouble( "slicedistance", 0 ) );
    
  this->Rewind();
  while ( this->Seek ( "image" ) ) 
    {
    tmpStr = this->ReadString( "name" );
    if ( tmpStr ) 
      {
      newStudy->push_back( tmpStr );
      }
    }
  
  if ( newStudy->GetMultiFile() )
    newStudy->SetDims( Width, Height, newStudy->size() );
  else
    newStudy->SetDims( Width, Height, Depth );
    
#ifdef DEBUG
  printf( "  Study describes %dx%dx%d data volume\n", Width, Height, newStudy->GetMultiFile() ? (int)newStudy->size() : newStudy->GetDims( AXIS_Z ) );
#endif
  
  if ( study ) delete study;
  study = newStudy;
  
  return *this;
}

Study*
ClassStream::ReadStudy( const char* dir )
{
  Study *study = NULL;

  ClassStream stream;
  stream.Open( dir, "images", ClassStream::READ );
  stream >> study;
  stream.Close();

  study->SetFileSystemPath( dir );
  return study;
}

} // namespace cmtk
