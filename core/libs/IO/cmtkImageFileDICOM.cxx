/*
//
//  Copyright 2004-2012 SRI International
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

#include <dcmtk/dcmimgle/diutils.h>

#ifndef DCM_RawDataType_ImageType
#define DCM_RawDataType_ImageType DcmTagKey(0x0043,0x102f)
#endif

#ifndef DCM_ManufacturerModelName
#define DCM_ManufacturerModelName DcmTagKey(0x0008,0x1090)
#endif

#ifndef DCM_PatientsName
#define DCM_PatientsName DCM_PatientName
#endif

namespace
cmtk
{

void
ImageFileDICOM::Print() const
{
  cmtk::DebugOutput( 1 ) << "  File Name =            [" << this->fpath << "/" << this->fname << "]\n";
  cmtk::DebugOutput( 1 ) << "  SeriesID =             [" << this->GetTagValue( DCM_SeriesInstanceUID ) << "]\n";
  cmtk::DebugOutput( 1 ) << "  StudyID =              [" << this->GetTagValue( DCM_StudyInstanceUID ) << "]\n";
  cmtk::DebugOutput( 1 ) << "  ImagePositionPatient = [" << this->GetTagValue( DCM_ImagePositionPatient ) << "]\n";
  cmtk::DebugOutput( 1 ) << "  AcquisitionNumber =    [" << this->AcquisitionNumber << "]\n";
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
  if ( this->IsMultislice || other.IsMultislice )
    return false;

  if ( ! disableCheckOrientation )
    {
    double orientThis[6], orientOther[6];
    sscanf( this->GetTagValue( DCM_ImageOrientationPatient ).c_str(), "%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf", orientThis, orientThis+1, orientThis+2, orientThis+3, orientThis+4, orientThis+5 );
    sscanf( other.GetTagValue( DCM_ImageOrientationPatient ).c_str(), "%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf", orientOther, orientOther+1, orientOther+2, orientOther+3, orientOther+4, orientOther+5 );

    for ( int i = 0; i < 6; ++i )
      {
      if ( fabs( orientThis[i] - orientOther[i] ) > numericalTolerance )
	return false;
      }
    }

  return
    ( this->GetTagValue( DCM_FrameOfReferenceUID ) == other.GetTagValue( DCM_FrameOfReferenceUID ) ) && 
    ( this->GetTagValue( DCM_SeriesInstanceUID ) == other.GetTagValue( DCM_SeriesInstanceUID ) ) && 
    ( this->GetTagValue( DCM_StudyInstanceUID ) == other.GetTagValue( DCM_StudyInstanceUID ) ) && 
    ( this->GetTagValue( DCM_EchoTime ) == other.GetTagValue( DCM_EchoTime ) ) &&
    ( this->GetTagValue( DCM_RepetitionTime ) == other.GetTagValue( DCM_RepetitionTime ) ) && 
    ( this->BValue == other.BValue ) &&
    ( this->BVector == other.BVector ) &&
    (( this->AcquisitionNumber == other.AcquisitionNumber ) || ignoreAcquisitionNumber) && 
    ( this->m_RawDataType == other.m_RawDataType );
}

ImageFileDICOM::ImageFileDICOM( const char* filename )
  :IsMultislice( false ),
   IsDWI( false ),   
   BValue( 0 ),
   BVector( 0.0 ),
   m_RawDataType( "unknown" )
{
  if ( cmtk::FileFormat::Identify( filename, false /*decompress*/ ) != cmtk::FILEFORMAT_DICOM ) // need to disable "decompress" in Identify() because DCMTK cannot currently read using on-the-fly decompression.
    throw(0);

  const char *last_slash = strrchr( filename, CMTK_PATH_SEPARATOR );
  if ( last_slash ) 
    {
    fname = strdup(last_slash+1);
    char *suffix = strrchr( fname, '.' );
    if ( suffix )
      if ( !strcmp( suffix, ".Z" ) || !strcmp( suffix, ".gz" ) ) 
	{
	*suffix = 0;
	}
    
    int path_len = last_slash-filename;
    fpath = (char*)malloc( path_len+1 );
    strncpy( fpath, filename, path_len );
    fpath[path_len] = 0;
    } 
  else
    {
    fname = strdup( filename );
    fpath = NULL;
    }
  
  std::auto_ptr<DcmFileFormat> fileformat( new DcmFileFormat );
  
  fileformat->transferInit();
  OFCondition status = fileformat->loadFile( filename );
  fileformat->transferEnd();
  
  if ( !status.good() ) 
    {
    cmtk::StdErr << "Error: cannot read DICOM file " << filename << " (" << status.text() << ")\n";
    throw (0);
    }
  
  DcmDataset *dataset = fileformat->getAndRemoveDataset();
  if ( ! dataset )
  {
     throw(1);
  }

  this->m_Document = std::auto_ptr<DiDocument>( new DiDocument( dataset, dataset->getOriginalXfer(), CIF_AcrNemaCompatibility ) );
  if ( ! this->m_Document.get() || ! this->m_Document->good() ) 
    {
    throw(2);
    }

  const char* tmpStr = NULL;

  if ( this->m_Document->getValue( DCM_Modality, tmpStr ) )
    this->m_TagToStringMap[DCM_Modality] = tmpStr;

  // check for multi-slice DICOMs
  Uint16 nFrames = 0;
  if ( this->m_Document->getValue( DCM_NumberOfFrames, nFrames ) ) 
    {
    this->IsMultislice = (nFrames > 1 );
    }

  if ( this->m_Document->getValue( DCM_PatientsName, tmpStr ) )
    this->m_TagToStringMap[DCM_PatientsName] = tmpStr;

  if ( this->m_Document->getValue( DCM_SeriesInstanceUID, tmpStr ) )
    this->m_TagToStringMap[DCM_SeriesInstanceUID] = tmpStr;

  if ( this->m_Document->getValue( DCM_FrameOfReferenceUID, tmpStr ) )
    this->m_TagToStringMap[DCM_FrameOfReferenceUID] = tmpStr;

  if ( this->m_Document->getValue( DCM_SeriesDescription, tmpStr ) )
    this->m_TagToStringMap[DCM_SeriesDescription] = tmpStr;

  if ( this->m_Document->getValue( DCM_StudyInstanceUID, tmpStr ) )
    this->m_TagToStringMap[DCM_StudyInstanceUID] = tmpStr;

  if ( this->m_Document->getValue( DCM_StudyID, tmpStr ) )
    this->m_TagToStringMap[DCM_StudyID] = tmpStr;

  if ( this->m_Document->getValue( DCM_StudyDate, tmpStr ) )
    this->m_TagToStringMap[DCM_StudyDate] = tmpStr;
  
  if ( this->m_Document->getValue( DCM_ImagePositionPatient, tmpStr ) )
    this->m_TagToStringMap[DCM_ImagePositionPatient] = tmpStr;

  if ( this->m_Document->getValue( DCM_ImageOrientationPatient, tmpStr ) )
    this->m_TagToStringMap[DCM_ImageOrientationPatient] = tmpStr;

  if ( ! this->m_Document->getValue( DCM_InstanceNumber, InstanceNumber ) )
    InstanceNumber = 0;

  if ( ! this->m_Document->getValue( DCM_AcquisitionNumber, AcquisitionNumber ) )
    AcquisitionNumber = 0;

  if ( this->m_Document->getValue( DCM_RescaleIntercept, tmpStr ) )
    this->m_TagToStringMap[DCM_RescaleIntercept] = tmpStr;

  if ( this->m_Document->getValue( DCM_RescaleSlope, tmpStr ) )
    this->m_TagToStringMap[DCM_RescaleSlope] = tmpStr;

  // check for MR modality and treat accordingly
  if ( this->GetTagValue( DCM_Modality ) == "MR" )
    {
    if ( this->m_Document->getValue( DCM_EchoTime, tmpStr ) )
      this->m_TagToStringMap[DCM_EchoTime] = tmpStr;
    
    if ( this->m_Document->getValue( DCM_RepetitionTime, tmpStr ) )
      this->m_TagToStringMap[DCM_RepetitionTime] = tmpStr;

    if ( this->m_Document->getValue( DCM_ImagingFrequency, tmpStr ) )
      this->m_TagToStringMap[DCM_ImagingFrequency] = tmpStr;
    }
  
  if ( this->m_Document->getValue( DCM_ManufacturerModelName, tmpStr ) != 0 )
    this->m_TagToStringMap[DCM_ManufacturerModelName] = tmpStr;

  if ( this->m_Document->getValue( DCM_DeviceSerialNumber, tmpStr ) != 0 )
    this->m_TagToStringMap[DCM_DeviceSerialNumber] = tmpStr;

  // check for which vendor and deal with specifics elsewhere
  if ( this->m_Document->getValue( DCM_Manufacturer, tmpStr ) != 0 )
    {
    this->m_TagToStringMap[DCM_Manufacturer] = tmpStr;
    if ( !strncmp( tmpStr, "SIEMENS", 7 ) )
      {
      this->DoVendorTagsSiemens();
      }      
    
    if ( !strncmp( tmpStr, "GE", 2 ) )
      {
      this->DoVendorTagsGE();
      }
    }
}

void
ImageFileDICOM::DoVendorTagsSiemens()
{
  Uint16 nFrames = 0;
  const char* tmpStr = NULL;

  this->IsMultislice = this->m_Document->getValue( DcmTagKey (0x0019,0x100a), nFrames ); // Number of Slices tag
  this->IsMultislice |= ( this->m_Document->getValue( DCM_ImageType, tmpStr ) && strstr( tmpStr, "MOSAIC" ) ); // mosaics are always multi-slice
  
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
    
    if ( (this->IsDWI = (this->m_Document->getValue( DcmTagKey(0x0019,0x100d), tmpStr )!=0)) ) // "Directionality" tag
      {
      if ( this->m_Document->getValue( DcmTagKey(0x0019,0x100c), tmpStr ) != 0 ) // bValue tag
	{
	this->BValue = atoi( tmpStr );
	this->IsDWI |= (this->BValue > 0);
	}
      
      if ( this->BValue > 0 )
	{
	for ( int idx = 0; idx < 3; ++idx )
	  {
	  this->IsDWI |= (this->m_Document->getValue( DcmTagKey(0x0019,0x100e), this->BVector[idx], idx ) != 0);
	  }
	}
      }
    }
}

void
ImageFileDICOM::DoVendorTagsGE()
{
  const char* tmpStr = NULL;
  int tmpInt = 0;

  if ( this->GetTagValue( DCM_Modality ) == "MR" )
    {
    // raw data type
    Sint16 rawTypeIdx = 3;
    if ( ! this->m_Document->getValue( DCM_RawDataType_ImageType, rawTypeIdx ) )
      rawTypeIdx = 0; // assume this is a magnitude image
    rawTypeIdx = std::min( 3, std::max( 0, (int)rawTypeIdx ) );
    
    const char *const RawDataTypeString[4] = { "magnitude", "phase", "real", "imaginary" };
    this->m_RawDataType = RawDataTypeString[rawTypeIdx];
    
    // dwi information
    this->IsDWI = false;
    if ( this->m_Document->getValue( DcmTagKey(0x0019,0x10e0), tmpStr ) > 0 ) // Number of Diffusion Directions
      {
      const int nDirections = atoi( tmpStr );
      if ( nDirections > 0 )
	{
	this->IsDWI = true;
	
	if ( this->m_Document->getValue( DcmTagKey(0x0043,0x1039), tmpStr ) > 0 ) // bValue tag
	  {
	  if ( 1 == sscanf( tmpStr, "%d\\%*c", &tmpInt ) )
	    {
	    this->BValue = static_cast<Sint16>( tmpInt );

	    for ( int i = 0; i < 3; ++i )
	      {
	      if ( this->m_Document->getValue( DcmTagKey(0x0019,0x10bb+i), tmpStr ) > 0 ) // bVector tags
		{
		this->BVector[i] = atof( tmpStr );
		}
	      else
		{
		this->BVector[i] = 0;
		}
	      }
	    this->BVector[2] *= -1; // for some reason z component apparently requires negation of sign on GE
	    }
	  }
	}
      }
    }
}

ImageFileDICOM::~ImageFileDICOM()
{
  free( fname );
  free( fpath );
}

bool
ImageFileDICOM::MatchAnyPattern( const std::map<DcmTagKey,std::string>& patterns ) const
{
  const char* tmpStr = NULL;

  // check for positive include list
  if ( !patterns.empty() )
    {
    for ( std::map<DcmTagKey,std::string>::const_iterator it = patterns.begin(); it != patterns.end(); ++it )
      {
      // if tag not found, do not include
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
  const char* tmpStr = NULL;

  // check for positive include list
  if ( !patterns.empty() )
    {
    for ( std::map<DcmTagKey,std::string>::const_iterator it = patterns.begin(); it != patterns.end(); ++it )
      {
      // if tag not found, do not include
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
