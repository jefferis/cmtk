/*
//
//  Copyright 1997-2012 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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
#include <System/cmtkDebugOutput.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkMetaInformationObject.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkStudy.h>
#include <IO/cmtkStudyImageSet.h>
#include <IO/cmtkVolumeFromStudy.h>
#include <IO/cmtkVolumeFromFile.h>
#include <IO/cmtkFileFormat.h>

#include <mxml.h>

#ifndef CMTK_USE_DCMTK
#error Build system is broken: this application should not be build if CMTK_USE_DCMTK is not set.
#endif

#ifdef _MSC_VER
#  include <windows.h>
#else
#  include <dirent.h>
#  include <fnmatch.h>
#  include <unistd.h>
#endif

#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <sstream>

#include <iostream>
#include <memory>
#include <fstream>

#include <dcmtk/dcmdata/dctk.h>
#include <dcmtk/dcmimgle/didocu.h>
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

const char* OutPathPattern = "%03d.nii";
std::vector<std::string> SearchRootDirVector;

std::ofstream cnull( "/dev/null" );

const char progress_chars[] = "-\\|/";
int progress = 0;

bool Recursive = false;
int SortFiles = 1;
bool WriteXML = false;

bool DisableOrientationCheck = false;
double Tolerance = 1e-5;

bool IgnoreAcquisitionNumber = false;

/// Enum type to select DICOM information to be embedded into output images as "description".
typedef enum
{
  /// No embedding.
  EMBED_NONE = 0,
  /// Embed StudyID plus Date.
  EMBED_STUDYID_STUDYDATE = 1,
  /// Embed patient name.
  EMBED_PATIENTNAME = 2,
  /// Embed series description.
  EMBED_SERIESDESCR = 3
} EmbedInfoEnum;

/// Selector for embedded image information.
EmbedInfoEnum EmbedInfo = EMBED_STUDYID_STUDYDATE;

/// Class handling a single image file and its meta data.
class ImageFile 
{
public:
  /// File name.
  char* fname;

  /// File system path (i.e., directory).
  char* fpath;

  /// Modality.
  std::string Modality;

  /// Device manufacturer.
  std::string Manufacturer;

  /// Device model.
  std::string ManufacturerModel;

  /// Flag for multislice images
  bool IsMultislice;

  /// Patient name.
  std::string PatientName;

  /// DICOM SeriesUID.
  std::string SeriesUID;

  /// DICOM SeriesDescription
  std::string SeriesDescription;

  /// DICOM StudyID.
  std::string StudyID;

  /// DICOM StudyDate.
  std::string StudyDate;
  
  /// DICOM StudyUID.
  std::string StudyUID;

  /// MR repetition time, TR
  std::string RepetitionTime;

  /// MR echo time, TE.
  std::string EchoTime;

  /// 3D image position (first pixel) in patient coordinates.
  std::string ImagePositionPatient;

  /// 3D image orientation (pos. x and pos. y image axes) in patient coordinates.
  std::string ImageOrientationPatient;

  /// DICOM acquisition number.
  Sint32 AcquisitionNumber;

  /// DICOM image number (index in volume).
  Sint32 InstanceNumber;

  /// Raw data type (real, imaginary, phase, magnitude) currently supported on GE images only.
  std::string RawDataType;

  /// Flag for diffusion-weighted images.
  bool IsDWI;

  /// B value for DWI.
  Sint16 BValue;

  /// B vector for DWI.
  cmtk::FixedVector<3,double> BVector;

  /// Constructor.
  ImageFile( const char* filename );

  /// Destructor.
  ~ImageFile();

  /// Determine whether two images match, i.e., belong to the same volume.
  bool Match( const ImageFile& other ) const;

  /// Compare order based on file name (for lexicographic sorting).
  static bool lessFileName( const ImageFile* lhs, const ImageFile* rhs )
  {
    return strcmp( lhs->fname, rhs->fname ) < 0;
  }

  /// Compare order based on image instace (for sorting in acquisition order).
  static bool lessInstanceNumber( const ImageFile* lhs, const ImageFile* rhs )
  {
    return lhs->InstanceNumber < rhs->InstanceNumber;
  }

  void Print() const;

private:
  /// Handle Siemens private tags.
  void DoVendorTagsSiemens( const DiDocument& document );

  /// Handle GE private tags.
  void DoVendorTagsGE( const DiDocument& document );
};

void
ImageFile::Print() const
{
  cmtk::DebugOutput( 1 ) << "  File Name =            [" << this->fpath << "/" << this->fname << "]\n";
  cmtk::DebugOutput( 1 ) << "  SeriesID =             [" << this->SeriesUID << "]\n";
  cmtk::DebugOutput( 1 ) << "  StudyID =              [" << this->StudyUID << "]\n";
  cmtk::DebugOutput( 1 ) << "  ImagePositionPatient = [" << this->ImagePositionPatient << "]\n";
  cmtk::DebugOutput( 1 ) << "  AcquisitionNumber =    [" << this->AcquisitionNumber << "]\n";
  cmtk::DebugOutput( 1 ) << "  Modality =             [" << this->Modality << "]\n";

  if ( this->Modality == "MR" )
    {
    cmtk::DebugOutput( 1 ) << "  EchoTime =          [" << this->EchoTime << "]\n";
    cmtk::DebugOutput( 1 ) << "  RepetitionTime =      [" << this->RepetitionTime << "]\n";
    }
}
  
bool
ImageFile::Match( const ImageFile& other ) const
{
  // do not stack multislice images
  if ( this->IsMultislice || other.IsMultislice )
    return false;

  if ( ! DisableOrientationCheck )
    {
    double orientThis[6], orientOther[6];
    sscanf( this->ImageOrientationPatient.c_str(), "%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf", orientThis, orientThis+1, orientThis+2, orientThis+3, orientThis+4, orientThis+5 );
    sscanf( other.ImageOrientationPatient.c_str(), "%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf", orientOther, orientOther+1, orientOther+2, orientOther+3, orientOther+4, orientOther+5 );

    for ( int i = 0; i < 6; ++i )
      {
      if ( fabs( orientThis[i] - orientOther[i] ) > Tolerance )
	return false;
      }
    }

  return
    ( SeriesUID == other.SeriesUID ) && 
    ( StudyUID == other.StudyUID ) && 
    ( EchoTime == other.EchoTime ) &&
    ( RepetitionTime == other.RepetitionTime ) && 
    ( BValue == other.BValue ) &&
    ( BVector == other.BVector ) &&
    (( AcquisitionNumber == other.AcquisitionNumber ) || IgnoreAcquisitionNumber) && 
    ( this->RawDataType == other.RawDataType );
}

ImageFile::ImageFile( const char* filename )
  : Modality( "unknown" ),
    Manufacturer( "unknown" ),
    ManufacturerModel( "unknown" ),
    IsMultislice( false ),
    RawDataType( "unknown" ),
    IsDWI( false ),
    BValue( 0 ),
    BVector( cmtk::FixedVector<3,double>::Init( 0.0 ) )
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

  std::auto_ptr<DiDocument> document( new DiDocument( dataset, dataset->getOriginalXfer(), CIF_AcrNemaCompatibility ) );
  if ( ! document.get() || ! document->good() ) 
    {
    throw(2);
    }

  const char* tmpStr = NULL;
  if ( document->getValue( DCM_Modality, tmpStr ) )
    this->Modality = tmpStr;

  // check for multi-slice DICOMs
  Uint16 nFrames = 0;
  if ( document->getValue( DCM_NumberOfFrames, nFrames ) ) 
    {
    this->IsMultislice = (nFrames > 1 );
    }

  if ( document->getValue( DCM_PatientsName, tmpStr ) )
    PatientName = tmpStr;

  if ( document->getValue( DCM_SeriesInstanceUID, tmpStr ) )
    SeriesUID = tmpStr;

  if ( document->getValue( DCM_SeriesDescription, tmpStr ) )
    SeriesDescription = tmpStr;

  if ( document->getValue( DCM_StudyInstanceUID, tmpStr ) )
    StudyUID = tmpStr;

  if ( document->getValue( DCM_StudyID, tmpStr ) )
    StudyID = tmpStr;

  if ( document->getValue( DCM_StudyDate, tmpStr ) )
    StudyDate = tmpStr;
  
  if ( document->getValue( DCM_ImagePositionPatient, tmpStr ) )
    ImagePositionPatient = tmpStr;

  if ( document->getValue( DCM_ImageOrientationPatient, tmpStr ) )
    ImageOrientationPatient = tmpStr;

  if ( ! document->getValue( DCM_InstanceNumber, InstanceNumber ) )
    InstanceNumber = 0;

  if ( ! document->getValue( DCM_AcquisitionNumber, AcquisitionNumber ) )
    AcquisitionNumber = 0;

  // check for MR modality and treat accordingly
  if ( this->Modality == "MR" )
    {
    if ( document->getValue( DCM_EchoTime, tmpStr ) )
      EchoTime = tmpStr;
    
    if ( document->getValue( DCM_RepetitionTime, tmpStr ) )
      RepetitionTime = tmpStr;
    }
  else
    {
    this->EchoTime = this->RepetitionTime = "";
    }
  
  if ( document->getValue( DCM_ManufacturerModelName, tmpStr ) != 0 )
    this->ManufacturerModel = tmpStr;

  // check for which vendor and deal with specifics elsewhere
  if ( document->getValue( DCM_Manufacturer, tmpStr ) != 0 )
    {
    this->Manufacturer = tmpStr;
    if ( !strncmp( tmpStr, "SIEMENS", 7 ) )
      {
      this->DoVendorTagsSiemens( *document );
      }      
    
    if ( !strncmp( tmpStr, "GE", 2 ) )
      {
      this->DoVendorTagsGE( *document );
      }
    }
}

void
ImageFile::DoVendorTagsSiemens( const DiDocument& document )
{
  Uint16 nFrames = 0;
  const char* tmpStr = NULL;

  this->IsMultislice = document.getValue( DcmTagKey (0x0019,0x100a), nFrames ); // Number of Slices tag
  
  if ( this->Modality == "MR" )
    {
    if ( (this->IsDWI = (document.getValue( DcmTagKey(0x0019,0x100d), tmpStr )!=0)) ) // "Directionality" tag
      {
      if ( document.getValue( DcmTagKey(0x0019,0x100c), tmpStr ) != 0 ) // bValue tag
	{
	this->BValue = atoi( tmpStr );
	this->IsDWI |= (this->BValue > 0);
	}
      
      if ( this->BValue > 0 )
	{
	for ( int idx = 0; idx < 3; ++idx )
	  {
	  this->IsDWI |= (document.getValue( DcmTagKey(0x0019,0x100e), this->BVector[idx], idx ) != 0);
	  }
	}
      }
    }
}

void
ImageFile::DoVendorTagsGE( const DiDocument& document )
{
  const char* tmpStr = NULL;
  int tmpInt = 0;

  if ( this->Modality == "MR" )
    {
    // raw data type
    Sint16 rawTypeIdx = 3;
    if ( ! document.getValue( DCM_RawDataType_ImageType, rawTypeIdx ) )
      rawTypeIdx = 3; // assume this is a magnitude image
    rawTypeIdx = std::min( 3, std::max( 0, (int)rawTypeIdx ) );
    
    const char *const RawDataTypeString[4] = { "magn", "phas", "real", "imag" };
    this->RawDataType = RawDataTypeString[rawTypeIdx];
    
    // dwi information
    this->IsDWI = false;
    if ( document.getValue( DcmTagKey(0x0019,0x10e0), tmpStr ) > 0 ) // Number of Diffusion Directions
      {
      const int nDirections = atoi( tmpStr );
      if ( nDirections > 0 )
	{
	this->IsDWI = true;
	
	if ( document.getValue( DcmTagKey(0x0043,0x1039), tmpStr ) > 0 ) // bValue tag
	  {
	  if ( 1 == sscanf( tmpStr, "%d\\%*c", &tmpInt ) )
	    {
	    this->BValue = static_cast<Sint16>( tmpInt );

	    for ( int i = 0; i < 3; ++i )
	      {
	      if ( document.getValue( DcmTagKey(0x0019,0x10bb+i), tmpStr ) > 0 ) // bVector tags
		{
		this->BVector[i] = atof( tmpStr );
		}
	      else
		{
		this->BVector[i] = 0;
		}
	      }
	    }
	  }
	}
      }
    }
}

ImageFile::~ImageFile()
{
  free( fname );
  free( fpath );
}

/// Class handling a stack of image files.
class ImageStack : public std::vector<ImageFile*> 
{
public:
  /// This class.
  typedef ImageStack Self;
  
  /// Add new DICOM image file to this stack.
  void AddImageFile( ImageFile *const image );

  /// Match new image file against this volume stack.
  bool Match ( const ImageFile *newImage ) const;

  /// Write XML sidecare file.
  void WriteXML ( const std::string& name ) const;

  /// Write to image file.
  void WriteImage ( const std::string& name ) const;

  /// Print stack information.
  void print() const;

private:
  /// Generate custom whitespaces for XML output.
  static const char *WhitespaceWriteMiniXML( mxml_node_t*, int where);

};

bool
ImageStack::Match ( const ImageFile *newImage ) const
{
  if ( empty() ) 
    return 1;

  const ImageFile *check = front();
  if ( check )
    {
    if ( !check->Match( *newImage ) )
      return 0;

    for ( const_iterator it = this->begin(); it != this->end(); ++it )
      {
      // if we already have an image in same location in this study, 
      // then bump to next study
      if ( (*it)->ImagePositionPatient == newImage->ImagePositionPatient )
	return 0;
      }
    return 1;
    }
  else
    return 0;
}

void
ImageStack::AddImageFile ( ImageFile *const newImage )
{
  iterator it = begin();
  for ( ; it != end(); ++it )
    if ( newImage->InstanceNumber < (*it)->InstanceNumber ) break;
  insert( it, newImage );
}

const char *
ImageStack::WhitespaceWriteMiniXML( mxml_node_t* node, int where)
{
  const char* name = node->value.element.name;
  
  typedef struct _wsLookupType
  {
    /// XML element name.
    const char* name;
    /// Table of whitespace sequences.
    const char* ws[4];
  } wsLookupType;

  static const wsLookupType wsLookup[] = 
  {
    { "manufacturer", { "\t", NULL, NULL, "\n" } },
    { "model", { "\t", NULL, NULL, "\n" } },
    { "tr", { "\t", NULL, NULL, "\n" } },
    { "te", { "\t", NULL, NULL, "\n" } },
    { "dwi", { "\t", "\n", "\t", "\n" } },
    { "bValue", { "\t\t", NULL, NULL, "\n" } },
    { "bVector", { "\t\t", NULL, NULL, "\n" } },
    { "dcmPath", { NULL, NULL, NULL, "\n" } },
    { "dcmFile", { "\t", NULL, NULL, "\n" } },
    { NULL, {NULL, NULL, NULL, NULL} }
  };

  if ( (where >= 0) && (where < 4) )
    {
    for ( size_t idx = 0; wsLookup[idx].name; ++idx )
      {
      if ( ! strcmp( name, wsLookup[idx].name ) )
	return wsLookup[idx].ws[where];
      }
    }

  switch ( where )
    {
    case MXML_WS_BEFORE_OPEN:
      return NULL;
    case MXML_WS_AFTER_OPEN:
      return "\n";
    case MXML_WS_BEFORE_CLOSE:
      return NULL;
    case MXML_WS_AFTER_CLOSE:
      return "\n";
    }

  return NULL;
}

// wrap tolower() - on Mac, system function is not compatible with std::transform()
static int cmtkWrapToLower( const int c )
{
  return tolower( c );
}

void
ImageStack::WriteXML( const std::string& fname ) const
{
  mxml_node_t *x_root = mxmlNewElement( NULL, "?xml version=\"1.0\" encoding=\"utf-8\"?" );

  mxml_node_t *x_device = mxmlNewElement( x_root, "device" );

  mxml_node_t *x_manufacturer = mxmlNewElement( x_device, "manufacturer" );
  mxmlNewText( x_manufacturer, 0, this->front()->Manufacturer.c_str() );
    
  mxml_node_t *x_model = mxmlNewElement( x_device, "model" );
  mxmlNewText( x_model, 0, this->front()->ManufacturerModel.c_str() );

  std::string modality = this->front()->Modality;
  std::transform( modality.begin(), modality.end(), modality.begin(), cmtkWrapToLower );
  
  mxml_node_t *x_modality = mxmlNewElement( x_root, modality.c_str() );
  if ( modality == "mr" )
    {
    mxml_node_t *x_tr = mxmlNewElement( x_modality, "tr");
    mxmlNewReal( x_tr, atof( this->front()->RepetitionTime.c_str() ) );
    
    mxml_node_t *x_te = mxmlNewElement( x_modality, "te");
    mxmlNewReal( x_te, atof( this->front()->EchoTime.c_str() ) );
    
    if ( this->front()->IsDWI )
      {
      mxml_node_t *x_dwi = mxmlNewElement( x_modality, "dwi" );
      
      mxml_node_t *x_bval = mxmlNewElement( x_dwi, "bValue");
      mxmlNewInteger( x_bval, this->front()->BValue );
      
      mxml_node_t *x_bvec = mxmlNewElement( x_dwi, "bVector");
      for ( int idx = 0; idx < 3; ++idx )
	{
	mxmlNewReal( x_bvec, this->front()->BVector[idx] );
	}
      }
    }
    
  mxml_node_t *x_dcmpath = mxmlNewElement( x_root, "dcmPath" );
  mxmlNewText( x_dcmpath, 0, this->front()->fpath );

  mxml_node_t *x_stack = mxmlNewElement( x_root, "stack" );

  for ( const_iterator it = this->begin(); it != this->end(); ++it ) 
    {
    mxml_node_t *x_dcmfile = mxmlNewElement( x_stack, "dcmFile" );
    mxmlNewText( x_dcmfile, 0, (*it)->fname );
    }

  FILE *file = fopen( fname.c_str(), "w" );
  if ( file )
    {
    mxmlSaveFile( x_root, file, Self::WhitespaceWriteMiniXML );
    fputs( "\n", file ); // end last line
    fclose( file );
    }
  else
    {
    cmtk::StdErr << "ERROR: could not open file " << fname << " for writing\n";
    }
  
  mxmlDelete( x_root );
}

void
ImageStack::WriteImage( const std::string& fname ) const
{
  const ImageFile *first = this->front();
    
  cmtk::UniformVolume::SmartPtr volume;
  if ( !first->IsMultislice )
    {
    cmtk::StudyImageSet studyImageSet;
    
    studyImageSet.SetImageFormat( cmtk::FILEFORMAT_DICOM );
    studyImageSet.SetImageDirectory( first->fpath );
    studyImageSet.SetMultiFile( true );
    
    for ( const_iterator it = this->begin(); it != this->end(); ++it ) 
      {
      studyImageSet.push_back( (*it)->fname );
      }
    
    volume = cmtk::VolumeFromStudy::Read( &studyImageSet );
    }
  else
    {
    char fullPath[PATH_MAX];
#ifdef MSC_VER
    snprintf( fullPath, sizeof( fullPath ), "%s\\%s", first->fpath, first->fname );
#else
    snprintf( fullPath, sizeof( fullPath ), "%s/%s", first->fpath, first->fname );
#endif

    volume = cmtk::VolumeFromFile::ReadDICOM( fullPath );
    }

  if ( volume )
    {
    switch ( EmbedInfo )
      {
      default:
      case EMBED_NONE:
	break;
      case EMBED_STUDYID_STUDYDATE:
	volume->SetMetaInfo( cmtk::META_IMAGE_DESCRIPTION, first->StudyID + "_" + first->StudyDate );
	break;
	break;
      case EMBED_PATIENTNAME:
	volume->SetMetaInfo( cmtk::META_IMAGE_DESCRIPTION, first->PatientName );
	break;
      case EMBED_SERIESDESCR:
	volume->SetMetaInfo( cmtk::META_IMAGE_DESCRIPTION, first->SeriesDescription );
	break;
      }
    
    cmtk::VolumeIO::Write( *volume, fname.c_str() );
    cmtk::DebugOutput( 1 ).GetStream().printf( "\nOutput file:%s\nImage size: %3dx%3dx%3d pixels\nPixel size: %.4fx%.4fx%.4f mm\n\n", 
					       fname.c_str(), volume->m_Dims[0], volume->m_Dims[1], volume->m_Dims[2], volume->m_Delta[0], volume->m_Delta[1], volume->m_Delta[2] );
    }
  else
    {
    cmtk::StdErr << "WARNING: No valid volume was read.\n";
    }
  
  cmtk::DebugOutput( 1 ) << "DICOM Information: \n"
			 << "  Description:   " << first->SeriesDescription << "\n"
			 << "  Series:        " << first->SeriesUID << "\n"
			 << "  Study:         " << first->StudyUID << "\n"
			 << "  Acquisition:   " << first->AcquisitionNumber << "\n"
			 << "  TR / TE:       " << first->RepetitionTime << "ms /" << first->EchoTime << "ms\n"
			 << "  Position:      " << first->ImagePositionPatient << "\n"
			 << "  Orientation:   " << first->ImageOrientationPatient << "\n"
			 << "  Raw Data Type: " << first->RawDataType << "\n";
    
  cmtk::DebugOutput( 1 ) << "\nImage List:\n";
  for ( const_iterator it = this->begin(); it != this->end(); ++it ) 
    {
    cmtk::DebugOutput( 1 ) << (*it)->fname << " ";
    }
  cmtk::DebugOutput( 1 ) << "\n====================================================\n";
}

/// Class handling a list of image stacks, i.e., volumes.
class VolumeList : 
  public std::vector<ImageStack*> 
{
public:
  /// Add a new image file to the correct stack or create a new stack.
  void AddImageFile( ImageFile *const newImage );
  
  /// Write all volumes in this list.
  void WriteVolumes();
};

inline std::string &
replacein(std::string &s, const std::string &sub, const std::string &other)
{
  assert(!sub.empty());
  for ( size_t b = s.find(sub, 0); b != s.npos; b = s.find(sub, b) )
    {
    s.replace(b, sub.size(), other);
    b += other.size();
    }
  return s;
}

/// Make a string legal in a path by replacing spaces and colons with "_".
std::string
MakeLegalInPath( const std::string& s )
{
  std::string result = s;

  result = replacein( result, " ", "_" );      
  result = replacein( result, ":", "_" );  

  return result;
}

void
VolumeList::WriteVolumes() 
{
  size_t cntSingleImages = 0;

  int idx = 1;
  std::map< std::string,std::vector<const ImageStack*> > pathToVolumeMap;
  for ( const_iterator it = begin(); it != end(); ++it ) 
    {
    if ( ((*it)->size() > 1) || (*(*it)->begin())->IsMultislice )
      {
      // replace place holders
      std::string path( OutPathPattern );
      replacein( path, "%D", MakeLegalInPath( (*it)[0][0]->SeriesDescription ) );
      replacein( path, "%R", MakeLegalInPath( (*it)[0][0]->RepetitionTime ) );
      replacein( path, "%E", MakeLegalInPath( (*it)[0][0]->EchoTime ) );
      replacein( path, "%T", (*it)[0][0]->RawDataType );
      
      if ( path.length() > PATH_MAX )
	cmtk::StdErr << "ERROR: output path exceeds maximum path length";
      else
	pathToVolumeMap[path].push_back( *it );
      }
    else
      {
      ++cntSingleImages;
      }
    }

  if ( cntSingleImages )
    {
    cmtk::StdErr << "WARNING: " << cntSingleImages << " single image(s) could not be assigned to multi-image stacks.\n";

    cmtk::DebugOutput( 1 ) << "\n====================================================\n";
    cmtk::DebugOutput( 1 ) << "WARNING: " << cntSingleImages << " single image(s) could not be assigned to multi-image stacks:\n\n";
    for ( const_iterator it = begin(); it != end(); ++it ) 
      {
      if ( ((*it)->size() == 1) && !(*(*it)->begin())->IsMultislice )
	{
	(*(*it))[0]->Print();
	cmtk::DebugOutput( 1 ) << "\n";
	}
      }
    cmtk::DebugOutput( 1 ) << "\n====================================================\n";
    }
  
  for ( std::map< std::string,std::vector<const ImageStack*> >::const_iterator it = pathToVolumeMap.begin(); it != pathToVolumeMap.end(); ++it )
    {						
    const size_t nVolumes = it->second.size();

    // if there is only one volume with the given output path, just write it
    if ( nVolumes == 1 )
      {					       
      // if there's a "number" tag, get rid of it.
      std::string uniquePath = it->first;
      replacein( uniquePath, "%n", "" );
      replacein( uniquePath, "%N", "" );

      char finalPath[PATH_MAX];
      sprintf( finalPath, uniquePath.c_str(), idx++ );
      it->second[0]->WriteImage( finalPath );

      if ( WriteXML )
	{
	strcat( finalPath, ".xml" );
	it->second[0]->WriteXML( finalPath );
	}
      }
    else
      {			
      // otherwise, make unique paths for each of them
      for ( size_t i = 0; i < nVolumes; ++i )
	{
	std::ostringstream numberString;
	numberString.width( 1 + static_cast<int>( log( (double)nVolumes ) / M_LN10 ) );
	numberString.fill( '0' );
	numberString << std::right << 1+i;

	std::string uniquePath = it->first;
	replacein( uniquePath, "%n", numberString.str() );
	replacein( uniquePath, "%N", "-" + numberString.str() );

	char finalPath[PATH_MAX];
	sprintf( finalPath, uniquePath.c_str(), idx++ );
	it->second[i]->WriteImage( finalPath );

	if ( WriteXML )
	  {
	  strcat( finalPath, ".xml" );
	  it->second[i]->WriteXML( finalPath );
	  }
	}
      }
    }
}


void
VolumeList::AddImageFile( ImageFile *const newImage )
{
  if ( empty() ) 
    {
    ImageStack *newImageStack = new ImageStack;
    newImageStack->AddImageFile( newImage );
    push_back( newImageStack );
    } 
  else
    {
    const_iterator it = begin();
    while ( it != end() ) 
      {
      ImageStack *study = *it;
      if ( study->Match( newImage ) ) 
	{
	study->AddImageFile( newImage );
	return;
	} 
      else 
	{
	++it;
	}
      }
    ImageStack *newImageStack = new ImageStack;
    newImageStack->AddImageFile( newImage );
    push_back( newImageStack );    
    }
}

int
traverse_directory( VolumeList& volumeList, const std::string& path, const char *wildcard )
{
  char fullname[PATH_MAX];

  std::vector<ImageFile*> fileList;

#ifdef _MSC_VER
  WIN32_FIND_DATA fData;
  char pattern[PATH_MAX];
  snprintf( pattern, sizeof( pattern ), "%s\\%s", path.c_str(), wildcard );
  HANDLE hFind = FindFirstFile( pattern, &fData);
  do
    {
    snprintf( fullname, sizeof( fullname ), "%s\\%s", path.c_str(), fData.cFileName );
    if ( fData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY )
      {
      if ( Recursive && (fData.cFileName[0] != '.') )
	{
	traverse_directory( volumeList, fullname, wildcard );
	}
      }
    else
      {
      try
	{
	fileList.push_back( new ImageFile( fullname ) );
	(cmtk::StdErr << "\r" << progress_chars[ ++progress % 4 ]).flush();
	}
      catch ( ... )
	{
	// not a valid DICOM file
	}
      }
    }
  while (FindNextFile(hFind, &fData) != 0);
#else    
  DIR *dir_pointer = opendir ( path.c_str() );
  if ( dir_pointer != NULL ) 
    {
    struct dirent *entry_pointer;

    while ( (entry_pointer = readdir(dir_pointer)) ) 
      {
      strcat( strcat( strcpy( fullname, path.c_str() ), "/"), entry_pointer->d_name );
      struct stat entry_status;
      if ( !stat(fullname, &entry_status) ) 
	{
	if ( S_ISDIR( entry_status.st_mode ) && Recursive && (entry_pointer->d_name[0] != '.') ) 
	  {
	  strcat( fullname, "/" );
	  traverse_directory( volumeList, fullname, wildcard );
	  } 
	else
	  {
	  if ( !fnmatch(wildcard,entry_pointer->d_name,FNM_PERIOD) ) 
	    {
	    try
	      {
	      fileList.push_back( new ImageFile( fullname ) );
	      (cmtk::StdErr << "\r" << progress_chars[ ++progress % 4 ]).flush();
	      }
	    catch ( ... )
	      {
	      // not a valid DICOM file
	      }
	    }
	  }
	}
      }
    (void) closedir(dir_pointer);
    }
#endif

  switch ( SortFiles )
    {
    case 0:
    default:
      break;
    case 1:
      std::sort( fileList.begin(), fileList.end(), ImageFile::lessFileName );
      break;
    case 2:
      std::sort( fileList.begin(), fileList.end(), ImageFile::lessInstanceNumber );
      break;
    }
  
  for ( std::vector<ImageFile*>::const_iterator it = fileList.begin(); it != fileList.end(); ++it )
    {
    try 
      {
      volumeList.AddImageFile( *it );
      }
    catch (int)
      {
      }
    }

  cmtk::StdErr << "\r";
  return 0;
}


int
doMain ( const int argc, const char *argv[] )
{
  if (! dcmDataDict.isDictionaryLoaded() ) 
    {
#ifdef CMTK_DCMDICTPATH
    if ( dcmDataDict.wrlock().loadDictionary( CMTK_DCMDICTPATH ) )
      {
      dcmDataDict.unlock();
      }
    else
#endif
#ifdef CMTK_DCMDICTPATH_INSTALL
    if ( dcmDataDict.wrlock().loadDictionary( CMTK_DCMDICTPATH_INSTALL ) )
      {
      dcmDataDict.unlock();
      }
    else
#endif
      {
      cmtk::StdErr << "Data dictionary not avaliable. Please set DCMDICTPATH variable as path to dicom.dic file.\n";
      throw cmtk::ExitException( 1 );
      }
    }

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "DICOM to Image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Combine sets of DICOM slices to 3D image stacks" );

    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Input", "Input Options");
    cl.AddSwitch( Key( 'r', "recurse" ), &Recursive, true, "Recurse into directories" );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options");
    cl.AddOption( Key( 'O', "out-pattern" ), &OutPathPattern, "Output image path pattern. Use the following substitutions: "
		  "printf-style %d variante (image number); "
		  "%n (image number with automatic number of digits); "
		  "%N (like %n, but with a hyphen '-' before number if there is more than one image); "
		  "%D (DICOM SeriesDescription); "
		  "%R (DICOM RepetitionTime - MRI only)"
		  "%E (DICOM EchoTime - MRI only)"
		  "%T (RawDataType - vendor-specific, currently GE MRI only)" );

    cl.AddSwitch( Key( 'x', "xml" ), &WriteXML, true, "Write XML sidecar file for each created image." );

    cmtk::CommandLine::EnumGroup<EmbedInfoEnum>::SmartPtr embedGroup = cl.AddEnum( "embed", &EmbedInfo, "Embed DICOM information into output images as 'description' (if supported by output file format)." );
    embedGroup->AddSwitch( Key( "StudyID_StudyDate" ), EMBED_STUDYID_STUDYDATE, "StudyID, tag (0020,0010), then underscore, followed by StudyDate, tag (0008,0020). "
			   "Date is appended because StudyID is four digits only and will repeat sooner or later." );
    embedGroup->AddSwitch( Key( "PatientName" ), EMBED_PATIENTNAME, "Patient name, tag (0010,0010)" );
    embedGroup->AddSwitch( Key( "SeriesDescription" ), EMBED_SERIESDESCR, "Series description, tag (0008,103e)" );
    embedGroup->AddSwitch( Key( "None" ), EMBED_NONE, "Embed no information - leave 'description' field empty." );
    cl.EndGroup();

    cl.BeginGroup( "Sorting", "Sorting Options")->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddSwitch( Key( "no-sort" ), &SortFiles, 0, "Do NOT sort files by file name (determines order when resolving spatial collisions)" );
    cl.AddSwitch( Key( "sort-name" ), &SortFiles, 1, "Sort files lexicographically by file name." );
    cl.AddSwitch( Key( "sort-instance" ), &SortFiles, 2, "Sort files by image instance number. Use this when file names are different lengths, etc." );
    cl.EndGroup();

    cl.BeginGroup( "Stacking", "Stacking Options")->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddSwitch( Key( "ignore-acq-number" ), &IgnoreAcquisitionNumber, true, "Ignore 'AcquisitionNumber' tag for image grouping, i.e., do not separate stacks based on this tag." );
    cl.AddSwitch( Key( "no-orientation-check" ), &DisableOrientationCheck, true, "Disable checking of image orientations (to avoid rounding issues)" );
    cl.AddOption( Key( "tolerance" ), &Tolerance, "Tolerance for floating-point comparisons (must be >= 0; 0 = exact matches only; default: 1e-5)" );
    cl.EndGroup();

    cl.AddParameterVector( &SearchRootDirVector, "SearchDirList", "List of directories to search for DICOM files. Subdirectories are also search if '--recurse' option is used." );
    
    cl.Parse( argc, const_cast<const char**>( argv ) );
    }
  catch ( const cmtk::CommandLine::Exception& ex )
    {
    cmtk::StdErr << ex << "\n";
    throw cmtk::ExitException( 1 );
    }
  
  if ( !dcmDataDict.rdlock().findEntry( "RawDataType_ImageType" ) )
    {
    dcmDataDict.unlock();
    dcmDataDict.wrlock().addEntry( new DcmDictEntry( 0x0043, 0x102f, EVR_SS, "RawDataType_ImageType", 1, 1, NULL, OFFalse, "GE" ) );
    dcmDataDict.unlock();
    }
  
  VolumeList volumeList;
  for ( std::vector<std::string>::const_iterator it = SearchRootDirVector.begin(); it != SearchRootDirVector.end(); ++it )
    {
    traverse_directory( volumeList, *it, "*" );
    }
  
  volumeList.WriteVolumes();

  return 0;
}

#include "cmtkSafeMain"
