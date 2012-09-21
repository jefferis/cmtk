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
#include <System/cmtkStrUtility.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkMetaInformationObject.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkStudy.h>
#include <IO/cmtkStudyImageSet.h>
#include <IO/cmtkVolumeFromStudy.h>
#include <IO/cmtkVolumeFromFile.h>
#include <IO/cmtkImageFileDICOM.h>

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
#include <fstream>

#include <dcmtk/dcmdata/dctk.h>

void
AddPatternToMap( std::map<DcmTagKey,std::string>& map, const char* pattern )
{
  const char* equals = strchr( pattern, '=' );
  if ( equals )
    {
    const std::string tagName( pattern, static_cast<size_t>( equals-pattern ) );

    const DcmDictEntry* dictEntry = dcmDataDict.rdlock().findEntry( tagName.c_str() );
    dcmDataDict.unlock();

    if ( dictEntry )
      {
      map.insert( std::pair<DcmTagKey,std::string>( dictEntry->getKey(), equals+1 ) );      
      }
    else
      {
      cmtk::StdErr << "WARNING: DCMTK data dictionary does not have a tag entry named '" << tagName << "' - ignoring this\n";
      }
    }
}

// Map of file includion patterns
std::map<DcmTagKey,std::string> IncludePatterns;

void
CallbackAddInclude( const char* pattern )
{
  AddPatternToMap( IncludePatterns, pattern );
}

// Map of file exclusion patterns
std::map<DcmTagKey,std::string> ExcludePatterns;

void
CallbackAddExclude( const char* pattern )
{
  AddPatternToMap( ExcludePatterns, pattern );
}

const char* OutPathPattern = "image%n.nii";
std::vector<std::string> SearchRootDirVector;

std::ofstream cnull( "/dev/null" );

const char progress_chars[] = "-\\|/";
int progress = 0;

bool Recursive = false;

/// Enum type to select sort order key.
typedef enum
{
  /// No sorting - take files in order stored in file system.
  SORT_NONE = 0,
  /// Sort by file name.
  SORT_FILENAME = 1,
  /// Sort by instance number.
  SORT_INSTANCE_NUMBER = 2
} SortKeyEnum;

SortKeyEnum SortFiles = SORT_INSTANCE_NUMBER;

bool WriteXML = false;

bool WriteSingleSlices = false;

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

using cmtk::ImageFileDICOM;

/// Class handling a stack of image files.
class ImageStack : public std::vector<ImageFileDICOM*> 
{
public:
  /// This class.
  typedef ImageStack Self;
  
  /// Add new DICOM image file to this stack.
  void AddImageFile( ImageFileDICOM *const image );

  /// Match new image file against this volume stack.
  bool Match ( const ImageFileDICOM *newImage ) const;

  /// Write XML sidecar file.
  void WriteXML ( const std::string& name /*!< Sidecar XML file name. */, const cmtk::UniformVolume& volume /*!< Previously written image volume - provides information about coordinate system etc. */ ) const;

  /// Write to image file.
  cmtk::UniformVolume::SmartConstPtr WriteImage ( const std::string& name ) const;

  /// Print stack information.
  void print() const;

private:
  /// Generate custom whitespaces for XML output.
  static const char *WhitespaceWriteMiniXML( mxml_node_t*, int where);

};

bool
ImageStack::Match ( const ImageFileDICOM *newImage ) const
{
  if ( empty() ) 
    return 1;

  const ImageFileDICOM *check = front();
  if ( check )
    {
    if ( !check->Match( *newImage, Tolerance, DisableOrientationCheck, IgnoreAcquisitionNumber ) )
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
ImageStack::AddImageFile ( ImageFileDICOM *const newImage )
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
    { "manufacturer",    { "\t", NULL, NULL, "\n" } },
    { "model",           { "\t", NULL, NULL, "\n" } },
    { "tr",              { "\t", NULL, NULL, "\n" } },
    { "te",              { "\t", NULL, NULL, "\n" } },
    { "type",            { "\t", NULL, NULL, "\n" } },
    { "dwi",             { "\t", "\n", "\t", "\n" } },
    { "bValue",          { "\t\t", NULL, NULL, "\n" } },
    { "bVector",         { "\t\t", NULL, NULL, "\n" } },
    { "bVectorImage",    { "\t\t", NULL, NULL, "\n" } },
    { "bVectorStandard", { "\t\t", NULL, NULL, "\n" } },
    { "dcmPath",         { NULL, NULL, NULL, "\n" } },
    { "dcmFile",         { "\t", NULL, NULL, "\n" } },
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
ImageStack::WriteXML( const std::string& fname, const cmtk::UniformVolume& volume ) const
{
  mxmlSetWrapMargin( 120 ); // make enough room for indented bVectorStandard
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

    if ( this->front()->RawDataType != "unknown" )
      {
      mxml_node_t *x_type = mxmlNewElement( x_modality, "type");
      mxmlNewText( x_type, 0, this->front()->RawDataType.c_str() );
      }
    
    if ( this->front()->IsDWI )
      {
      mxml_node_t *x_dwi = mxmlNewElement( x_modality, "dwi" );
      
      mxml_node_t *x_bval = mxmlNewElement( x_dwi, "bValue");
      mxmlNewInteger( x_bval, this->front()->BValue );
      
      mxml_node_t *x_bvec = mxmlNewElement( x_dwi, "bVector");
      mxmlElementSetAttr( x_bvec, "coordinateSpace", "LPS" );
      for ( size_t idx = 0; idx < 3; ++idx )
	{
	mxmlNewReal( x_bvec, this->front()->BVector[idx] );
	}

      // Determine bVector in image LPS coordinate space:
      // First, create copy of image grid
      cmtk::UniformVolume::SmartPtr gridLPS = volume.CloneGrid();
      // Make sure still in LPS DICOM coordinate space
      gridLPS->ChangeCoordinateSpace( "LPS" );
      // Apply inverse of remaining image-to-space matrix to original bVector
      const cmtk::UniformVolume::CoordinateVectorType bVectorImage = this->front()->BVector * cmtk::Matrix3x3<cmtk::Types::Coordinate>( gridLPS->GetImageToPhysicalMatrix().GetInverse() );
      
      mxml_node_t *x_bvec_image = mxmlNewElement( x_dwi, "bVectorImage");
      mxmlElementSetAttr( x_bvec_image, "imageOrientation", gridLPS->GetMetaInfo( cmtk::META_IMAGE_ORIENTATION ).c_str() );
      for ( size_t idx = 0; idx < 3; ++idx )
	{
	mxmlNewReal( x_bvec_image, bVectorImage[idx] );
	}

      // Determine bVector in image RAS standard coordinate space:
      // First, create copy of image grid
      cmtk::UniformVolume::SmartPtr gridRAS = gridLPS->GetReoriented();
      // Apply inverse of remaining image-to-space matrix to original bVector
      const cmtk::UniformVolume::CoordinateVectorType bVectorStandard = this->front()->BVector * cmtk::Matrix3x3<cmtk::Types::Coordinate>( gridRAS->GetImageToPhysicalMatrix().GetInverse() );
      
      mxml_node_t *x_bvec_std = mxmlNewElement( x_dwi, "bVectorStandard");
      mxmlElementSetAttr( x_bvec_std, "imageOrientation", gridRAS->GetMetaInfo( cmtk::META_IMAGE_ORIENTATION ).c_str() );
      for ( size_t idx = 0; idx < 3; ++idx )
	{
	mxmlNewReal( x_bvec_std, bVectorStandard[idx] );
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

cmtk::UniformVolume::SmartConstPtr
ImageStack::WriteImage( const std::string& fname ) const
{
  const ImageFileDICOM *first = this->front();
    
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

  return volume;
}

/// Class handling a list of image stacks, i.e., volumes.
class VolumeList : 
  public std::vector<ImageStack*> 
{
public:
  /// Add a new image file to the correct stack or create a new stack.
  void AddImageFile( ImageFileDICOM *const newImage );
  
  /// Write all volumes in this list.
  void WriteVolumes();
};

std::string
PutNumberAndSanitize( const std::string& path, const std::string& numberString )
{
  return cmtk::StrReplace( cmtk::StrReplace( path, "%n", numberString ), "%N", "-" + numberString );
}

void
VolumeList::WriteVolumes() 
{
  size_t cntSingleImages = 0;

  std::map< std::string,std::vector<const ImageStack*> > pathToVolumeMap;
  for ( const_iterator it = begin(); it != end(); ++it ) 
    {
    if ( ((*it)->size() > 1) || (*(*it)->begin())->IsMultislice || WriteSingleSlices )
      {
      // replace place holders
      std::string path( OutPathPattern );
      path = cmtk::StrReplace( path, "%D", cmtk::StrMakeLegalInPath( (*it)[0][0]->SeriesDescription ) );
      path = cmtk::StrReplace( path, "%R", cmtk::StrMakeLegalInPath( (*it)[0][0]->RepetitionTime ) );
      path = cmtk::StrReplace( path, "%E", cmtk::StrMakeLegalInPath( (*it)[0][0]->EchoTime ) );
      path = cmtk::StrReplace( path, "%T", (*it)[0][0]->RawDataType );
      
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
      std::string uniquePath = PutNumberAndSanitize( it->first, "" );

      cmtk::UniformVolume::SmartConstPtr volume = it->second[0]->WriteImage( uniquePath.c_str() );

      if ( WriteXML )
	{
	it->second[0]->WriteXML( (uniquePath+".xml").c_str(), *volume );
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

	std::string uniquePath = PutNumberAndSanitize( it->first, numberString.str() );

	cmtk::UniformVolume::SmartConstPtr volume = it->second[i]->WriteImage( uniquePath.c_str() );

	if ( WriteXML )
	  {
	  it->second[i]->WriteXML( (uniquePath + ".xml").c_str(), *volume );
	  }
	}
      }
    }
}


void
VolumeList::AddImageFile( ImageFileDICOM *const newImage )
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

  std::vector<ImageFileDICOM*> fileList;

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
	ImageFileDICOM* newFile = new ImageFileDICOM( fullname );

	if ( newFile->MatchAllPatterns( IncludePatterns ) && !newFile->MatchAnyPattern( ExcludePatterns ) )
	  {
	  newFile->ReleaseDocument();
	  fileList.push_back( newFile );
	  }
	else
	  {
	  delete newFile;
	  }

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
	      ImageFileDICOM* newFile = new ImageFileDICOM( fullname );
	      
	      if ( newFile->MatchAllPatterns( IncludePatterns ) && !newFile->MatchAnyPattern( ExcludePatterns ) )
		{
		newFile->ReleaseDocument();
		fileList.push_back( new ImageFileDICOM( fullname ) );
		}
	      else
		{
		delete newFile;
		}
	      
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
    case SORT_NONE:
    default:
      break;
    case SORT_FILENAME:
      std::sort( fileList.begin(), fileList.end(), ImageFileDICOM::lessFileName );
      break;
    case SORT_INSTANCE_NUMBER:
      std::sort( fileList.begin(), fileList.end(), ImageFileDICOM::lessInstanceNumber );
      break;
    }
  
  for ( std::vector<ImageFileDICOM*>::const_iterator it = fileList.begin(); it != fileList.end(); ++it )
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
		  "%R (DICOM RepetitionTime - MRI only); "
		  "%E (DICOM EchoTime - MRI only); "
		  "%T (RawDataType - vendor-specific, currently GE MRI only)" );

    cl.AddSwitch( Key( 'x', "xml" ), &WriteXML, true, "Write XML sidecar file for each created image." );

    cmtk::CommandLine::EnumGroup<EmbedInfoEnum>::SmartPtr embedGroup = cl.AddEnum( "embed", &EmbedInfo, "Embed DICOM information into output images as 'description' (if supported by output file format)." );
    embedGroup->AddSwitch( Key( "StudyID_StudyDate" ), EMBED_STUDYID_STUDYDATE, "StudyID, tag (0020,0010), then underscore, followed by StudyDate, tag (0008,0020). "
			   "Date is appended because StudyID is four digits only and will repeat sooner or later." );
    embedGroup->AddSwitch( Key( "PatientName" ), EMBED_PATIENTNAME, "Patient name, tag (0010,0010)" );
    embedGroup->AddSwitch( Key( "SeriesDescription" ), EMBED_SERIESDESCR, "Series description, tag (0008,103e)" );
    embedGroup->AddSwitch( Key( "None" ), EMBED_NONE, "Embed no information - leave 'description' field empty." );
    cl.EndGroup();

    cl.BeginGroup( "Filtering", "Filtering Options")->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddCallback( Key( "filter" ), &CallbackAddInclude, "Filter DICOM files and include only those matching the given pattern of the form 'TagName=text', such that the value of the DICOM tag with the given name contains the given "
		    "text. If multiple filter patterns are provided via repeated use of this option, only files that match ALL patterns are included." );
    cl.AddCallback( Key( "exclude" ), &CallbackAddExclude, "Exclude all DICOM files matching the given pattern of the form 'TagName=text', such that the value of the DICOM tag with the given name contains the given text. "
      "If multiple exclusion patterns are provided, all files are excluded that match ANY of the patterns." );
    cl.EndGroup();

    cl.BeginGroup( "Sorting", "Sorting Options")->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddSwitch( Key( "no-sort" ), &SortFiles, SORT_NONE, "Do NOT sort files by file name (sorting determines image stack order when resolving spatial collisions)" );
    cl.AddSwitch( Key( "sort-by-name" ), &SortFiles, SORT_FILENAME, "Sort files lexicographically by file name. Use this when instance numbers are non-unique." );
    cl.AddSwitch( Key( "sort-by-instance" ), &SortFiles, SORT_INSTANCE_NUMBER, "Sort files by image instance number. Use this when file names are different lengths, etc." );
    cl.EndGroup();

    cl.BeginGroup( "Stacking", "Stacking Options")->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddSwitch( Key( "write-single-slices" ), &WriteSingleSlices, true, "Also write output images for single-slice DICOM files that could not be assigned to any 3D stacks. By default, these are skipped." );
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
