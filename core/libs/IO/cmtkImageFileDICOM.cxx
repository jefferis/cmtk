/*
//
//  Copyright 2004-2013 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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
//  $Revision: 4497 $
//
//  $LastChangedDate: 2012-08-24 13:46:21 -0700 (Fri, 24 Aug 2012) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include "cmtkImageFileDICOM.h"

#include <System/cmtkDebugOutput.h>

#include <IO/cmtkFileFormat.h>
#include <IO/cmtkSiemensCSAHeader.h>

#include <dcmtk/dcmimgle/diutils.h>

#include <sstream>

namespace
cmtk
{

void
ImageFileDICOM::Print() const
{
  cmtk::DebugOutput( 1 ) << "  File Name =            [" << this->m_FileDir << "/" << this->m_FileName << "]\n";
  cmtk::DebugOutput( 1 ) << "  SeriesID =             [" << this->GetTagValue( DCM_SeriesInstanceUID ) << "]\n";
  cmtk::DebugOutput( 1 ) << "  StudyID =              [" << this->GetTagValue( DCM_StudyInstanceUID ) << "]\n";
  cmtk::DebugOutput( 1 ) << "  ImagePositionPatient = [" << this->GetTagValue( DCM_ImagePositionPatient ) << "]\n";
  cmtk::DebugOutput( 1 ) << "  AcquisitionNumber =    [" << this->m_AcquisitionNumber << "]\n";
  cmtk::DebugOutput( 1 ) << "  Modality =             [" << this->GetTagValue( DCM_Modality ) << "]\n";

  if ( this->GetTagValue( DCM_Modality ) == "MR" )
    {
    cmtk::DebugOutput( 1 ) << "  EchoTime =          [" << this->GetTagValue( DCM_EchoTime ) << "]\n";
    cmtk::DebugOutput( 1 ) << "  RepetitionTime =      [" << this->GetTagValue( DCM_RepetitionTime ) << "]\n";
    }
}
  
bool
ImageFileDICOM::Match( const Self& other, const Types::Coordinate numericalTolerance, const bool disableCheckOrientation, const bool ignoreAcquisitionNumber ) const
{
  // do not stack multislice images
  if ( this->m_IsMultislice || other.m_IsMultislice )
    return false;

  if ( ! disableCheckOrientation )
    {
    double orientThis[6];
    if ( 6 != sscanf( this->GetTagValue( DCM_ImageOrientationPatient ).c_str(), "%30lf%*c%30lf%*c%30lf%*c%30lf%*c%30lf%*c%30lf", orientThis, orientThis+1, orientThis+2, orientThis+3, orientThis+4, orientThis+5 ) )
      return false;

    double orientOther[6];
    if ( 6 != sscanf( other.GetTagValue( DCM_ImageOrientationPatient ).c_str(), "%30lf%*c%30lf%*c%30lf%*c%30lf%*c%30lf%*c%30lf", orientOther, orientOther+1, orientOther+2, orientOther+3, orientOther+4, orientOther+5 ) )
      return false;

    for ( int i = 0; i < 6; ++i )
      {
      if ( fabs( orientThis[i] - orientOther[i] ) > numericalTolerance )
	return false;
      }
    }

  return
    ( this->m_FileDir == other.m_FileDir ) && // do not stack across directories
    ( this->GetTagValue( DCM_FrameOfReferenceUID ) == other.GetTagValue( DCM_FrameOfReferenceUID ) ) && 
    ( this->GetTagValue( DCM_SeriesInstanceUID ) == other.GetTagValue( DCM_SeriesInstanceUID ) ) && 
    ( this->GetTagValue( DCM_StudyInstanceUID ) == other.GetTagValue( DCM_StudyInstanceUID ) ) && 
    ( this->GetTagValue( DCM_EchoTime ) == other.GetTagValue( DCM_EchoTime ) ) &&
    ( this->GetTagValue( DCM_RepetitionTime ) == other.GetTagValue( DCM_RepetitionTime ) ) && 
    ( this->GetTagValue( DCM_InversionTime ) == other.GetTagValue( DCM_InversionTime ) ) && 
    ( this->m_BValue == other.m_BValue ) &&
    ( this->m_BVector == other.m_BVector ) &&
    (( this->m_AcquisitionNumber == other.m_AcquisitionNumber ) || ignoreAcquisitionNumber) && 
    ( this->m_RawDataType == other.m_RawDataType );
}

ImageFileDICOM::ImageFileDICOM( const std::string& filepath )
  :m_IsMultislice( false ),
   m_IsDWI( false ),   
   m_BValue( 0 ),
   m_BVector( 0.0 ),
   m_RawDataType( "unknown" )
{
  if ( cmtk::FileFormat::Identify( filepath, false /*decompress*/ ) != cmtk::FILEFORMAT_DICOM ) // need to disable "decompress" in Identify() because DCMTK cannot currently read using on-the-fly decompression.
    throw(0);

  this->m_FileName = filepath;
  this->m_FileDir = "";

  const size_t lastSlash = this->m_FileName.rfind( CMTK_PATH_SEPARATOR );
  if ( lastSlash != std::string::npos ) 
    {
    this->m_FileDir = this->m_FileName.substr( 0, lastSlash );

    // trim trailing slashes, both forward and back
    const size_t lastNotSlash = this->m_FileDir.find_last_not_of( "/\\" );
    if ( lastNotSlash != std::string::npos )
      {
      this->m_FileDir.erase( lastNotSlash+1 );
      }
    else
      {
      this->m_FileDir.clear();
      }

    this->m_FileName = this->m_FileName.substr( lastSlash+1 );

    const size_t suffix = this->m_FileName.rfind( '.' );
    if ( suffix != std::string::npos )
      {
      const std::string suffixStr = this->m_FileName.substr( suffix+1 );
      if ( (suffixStr == ".Z") || (suffixStr == ".gz") ) 
	{
	this->m_FileName = this->m_FileName.erase( suffix );
	}
      }

    } 
  
  std::auto_ptr<DcmFileFormat> fileformat( new DcmFileFormat );
  
  fileformat->transferInit();
  OFCondition status = fileformat->loadFile( filepath.c_str() );
  fileformat->transferEnd();
  
  if ( !status.good() ) 
    {
    cmtk::StdErr << "Error: cannot read DICOM file " << filepath << " (" << status.text() << ")\n";
    throw (0);
    }
  
  this->m_Dataset = std::auto_ptr<DcmDataset>( fileformat->getAndRemoveDataset() );
  if ( !this->m_Dataset.get() )
    {
    throw(1);
    }
  
  this->m_Document = std::auto_ptr<DiDocument>( new DiDocument( this->m_Dataset.get(), this->m_Dataset->getOriginalXfer(), CIF_AcrNemaCompatibility ) );
  if ( ! this->m_Document.get() || ! this->m_Document->good() ) 
    {
    throw(2);
    }

  // read the string tags that we need for later
  const DcmTagKey defaultStringTags[] = { DCM_Manufacturer, DCM_ManufacturerModelName, DCM_DeviceSerialNumber, DCM_StationName,
					  DCM_PatientsName, 
					  DCM_Modality, DCM_EchoTime, DCM_RepetitionTime, DCM_InversionTime, DCM_ImagingFrequency, DCM_SequenceName,
					  DCM_GE_PulseSequenceName, DCM_GE_PulseSequenceDate, DCM_GE_InternalPulseSequenceName,
					  DCM_StudyInstanceUID, DCM_StudyID, DCM_StudyDate,
					  DCM_FrameOfReferenceUID, DCM_SeriesInstanceUID, DCM_SeriesDescription,
					  DCM_ImagePositionPatient, DCM_ImageOrientationPatient, 
					  DCM_RescaleIntercept, DCM_RescaleSlope,
					  DcmTagKey( 0, 0 ) };
  
  for ( size_t tagIdx = 0; (defaultStringTags[tagIdx].getGroup() != 0) && (defaultStringTags[tagIdx].getElement() != 0); ++tagIdx )
    {
    const char* tmpStr = NULL;
    if ( this->m_Document->getValue( defaultStringTags[tagIdx], tmpStr ) )
      this->m_TagToStringMap[defaultStringTags[tagIdx]] = tmpStr;
    }

  // check for multi-slice DICOMs
  Uint16 nFrames = 0;
  if ( this->m_Document->getValue( DCM_NumberOfFrames, nFrames ) ) 
    {
    this->m_IsMultislice = (nFrames > 1 );
    }

  // read integer tags that we also need
  if ( ! this->m_Document->getValue( DCM_InstanceNumber, this->m_InstanceNumber ) )
    this->m_InstanceNumber = 0;

  if ( ! this->m_Document->getValue( DCM_AcquisitionNumber, this->m_AcquisitionNumber ) )
    this->m_AcquisitionNumber = 0;

  // check for which vendor and deal with specifics elsewhere
  if ( this->m_TagToStringMap[DCM_Manufacturer].substr( 0, 7 ) == "SIEMENS" )
    {
    this->DoVendorTagsSiemens();
    }      
  
  if ( this->m_TagToStringMap[DCM_Manufacturer].substr( 0, 2 ) == "GE" )
    {
    this->DoVendorTagsGE();
    }
}

void
ImageFileDICOM::DoVendorTagsSiemens()
{
  Uint16 nFrames = 0;
  const char* tmpStr = NULL;

  this->m_IsMultislice = (0 != this->m_Document->getValue( DcmTagKey (0x0019,0x100a), nFrames )); // Number of Slices tag
  this->m_IsMultislice |= ( this->m_Document->getValue( DCM_ImageType, tmpStr ) && strstr( tmpStr, "MOSAIC" ) ); // mosaics are always multi-slice
  
  if ( this->GetTagValue( DCM_Modality ) == "MR" )
    {
    // try to extract raw data type from "ImageType"
    if ( this->m_Document->getValue( DCM_ImageType, tmpStr ) )
      {
      if ( strstr( tmpStr, "\\P\\" ) )
	this->m_RawDataType = "phase";
      else if ( strstr( tmpStr, "\\M\\" ) )
	this->m_RawDataType = "magnitude";
      else if ( strstr( tmpStr, "\\R\\" ) )
	this->m_RawDataType = "real";
      }

    // for DWI, first check standard DICOM vendor tags
    if ( (this->m_IsDWI = (this->m_Document->getValue( DcmTagKey(0x0019,0x100d), tmpStr )!=0)) ) // "Directionality" tag
      {
      if ( this->m_Document->getValue( DcmTagKey(0x0019,0x100c), tmpStr ) != 0 ) // bValue tag
	{
	this->m_BValue = atoi( tmpStr );
	this->m_IsDWI |= (this->m_BValue > 0);
	}
      
      if ( this->m_BValue > 0 )
	{
	for ( int idx = 0; idx < 3; ++idx )
	  {
	  this->m_IsDWI |= (this->m_Document->getValue( DcmTagKey(0x0019,0x100e), this->m_BVector[idx], idx ) != 0);
	  }
	}
      }
    else
      {
      // no hint of DWI in standard tags, look into CSA header to confirm
      const Uint8* csaHeaderInfo = NULL;
      unsigned long csaHeaderLength = 0;
      if ( this->m_Dataset->findAndGetUint8Array ( DcmTagKey(0x0029,0x1010), csaHeaderInfo, &csaHeaderLength ).status() == OF_ok ) // this is expected in the "Image" CSA header, not the "Series" header.
	{
	SiemensCSAHeader csaHeader( (const char*)csaHeaderInfo, csaHeaderLength );
	SiemensCSAHeader::const_iterator it = csaHeader.find( "DiffusionDirectionality" );
	if ( (it != csaHeader.end()) && !it->second.empty() )
	  {
	  this->m_IsDWI = (0 == it->second[0].compare( 0, 11, "DIRECTIONAL" ));
	  }

	it = csaHeader.find( "B_value" );
	if ( (it != csaHeader.end()) && !it->second.empty() )
	  {
	  this->m_BValue = atof( it->second[0].c_str() );
	  }

	it = csaHeader.find( "DiffusionGradientDirection" );	
	if ( (it != csaHeader.end()) && (it->second.size() >= 3) )
	  {
	  for ( int idx = 0; idx < 3; ++idx )
	    {
	    this->m_BVector[idx] = atof( it->second[idx].c_str() );
	    }
	  }
	}
      }
    }
}

void
ImageFileDICOM::DoVendorTagsGE()
{
  int tmpInt = 0;

  if ( this->GetTagValue( DCM_Modality ) == "MR" )
    {
    // raw data type
    Sint16 rawTypeIdx = 3;
    if ( ! this->m_Document->getValue( DCM_GE_RawDataType_ImageType, rawTypeIdx ) )
      rawTypeIdx = 0; // assume this is a magnitude image
    rawTypeIdx = std::min( 3, std::max( 0, (int)rawTypeIdx ) );
    
    const char *const RawDataTypeString[4] = { "magnitude", "phase", "real", "imaginary" };
    this->m_RawDataType = RawDataTypeString[rawTypeIdx];

    Sint16 effEchoSpacing = 0;
    if ( this->m_Document->getValue( DCM_GE_EffectiveEchoSpacing, effEchoSpacing ) )
      {
      std::ostringstream ss;
      ss << effEchoSpacing;
      this->m_TagToStringMap[DCM_GE_EffectiveEchoSpacing] = ss.str();
      }
    
    // dwi information
    this->m_IsDWI = false;
    const char* tmpStr = NULL;
    if ( this->m_Document->getValue( DcmTagKey(0x0019,0x10e0), tmpStr ) > 0 ) // Number of Diffusion Directions
      {
      const int nDirections = atoi( tmpStr );
      if ( nDirections > 0 )
	{
	this->m_IsDWI = true;
	
	if ( this->m_Document->getValue( DcmTagKey(0x0043,0x1039), tmpStr ) > 0 ) // bValue tag
	  {
	  if ( 1 == sscanf( tmpStr, "%6d\\%*c", &tmpInt ) )
	    {
	    this->m_BValue = static_cast<Sint16>( tmpInt );

	    for ( int i = 0; i < 3; ++i )
	      {
	      if ( this->m_Document->getValue( DcmTagKey(0x0019,0x10bb+i), tmpStr ) > 0 ) // bVector tags
		{
		this->m_BVector[i] = atof( tmpStr );
		}
	      else
		{
		this->m_BVector[i] = 0;
		}
	      }
	    this->m_BVector[2] *= -1; // for some reason z component apparently requires negation of sign on GE
	    }
	  }
	}
      }
    }
}

bool
ImageFileDICOM::MatchAnyPattern( const std::map<DcmTagKey,std::string>& patterns ) const
{
  // check for positive include list
  if ( !patterns.empty() )
    {
    for ( std::map<DcmTagKey,std::string>::const_iterator it = patterns.begin(); it != patterns.end(); ++it )
      {
      // if tag not found, do not include
      const char* tmpStr = NULL;
      if ( this->m_Document->getValue( it->first, tmpStr ) )
	{
	// if tag value matches, then return true
	if ( strstr( tmpStr, it->second.c_str() ) )
	  return true;
	}
      }
    }

  return false; // did not find a match or list was empty
}

bool
ImageFileDICOM::MatchAllPatterns( const std::map<DcmTagKey,std::string>& patterns ) const
{
  // check for positive include list
  if ( !patterns.empty() )
    {
    for ( std::map<DcmTagKey,std::string>::const_iterator it = patterns.begin(); it != patterns.end(); ++it )
      {
      // if tag not found, do not include
      const char* tmpStr = NULL;
      if ( this->m_Document->getValue( it->first, tmpStr ) )
	{
	// if tag value matches, then return true
	if ( !strstr( tmpStr, it->second.c_str() ) )
	  return false;
	}
      }
    }

  return true; // all matched or list empty
}


} // namespace CMTK
