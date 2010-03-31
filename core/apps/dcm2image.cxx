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

#include <cmtkconfig.h>

#include <cmtkConsole.h>
#include <cmtkCommandLine.h>

#include <cmtkVolumeIO.h>
#include <cmtkStudy.h>
#include <cmtkStudyImageSet.h>
#include <cmtkVolumeFromStudy.h>
#include <cmtkFileFormat.h>

#ifndef CMTK_HAVE_DCMTK
#error Build system is broken: this application should not be build if CMTK_HAVE_DCMTK is not set.
#endif

#ifdef _MSC_VER
#  include <windows.h>
#else
#  include <dirent.h>
#  include <fnmatch.h>
#  include <unistd.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>

#include <vector>
#include <map>
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

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace dcm2image
{
#endif
const char* OutPathPattern = "%03d.hdr";
const char* SearchRootDir = ".";

std::ofstream cnull( "/dev/null" );

const char progress_chars[] = "-\\|/";
int progress = 0;

bool Recursive = false;
bool SortFilesByName = true;

bool Verbose = false;

bool WithExtensionsGE = false;
const char *const GERawDataTypeString[4] = { "magn", "phas", "real", "imag" };

bool DisableOrientationCheck = false;
double Tolerance = 1e-5;

class ImageFileDCM 
{
public:
  /// File name.
  char* fname;

  /// File system path (i.e., directory).
  char* fpath;
  
  /// DICOM SeriesUID.
  std::string SeriesUID;

  /// DICOM SeriesDescription
  std::string SeriesDescription;

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
  Sint32 ImageNumber;

  /// GE private DICOM tag: raw data type (real, imaginary, phase, magnitude).
  Sint16 GERawDataType;

  /// Constructor.
  ImageFileDCM( const char* filename );

  /// Destructor.
  ~ImageFileDCM();

  /// Determine whether two images match, i.e., belong to the same volume.
  bool Match( const ImageFileDCM& other ) const;
};

bool
ImageFileDCM::Match( const ImageFileDCM& other ) const
{
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
    ( SeriesUID == other.SeriesUID ) && ( StudyUID == other.StudyUID ) && 
    ( EchoTime == other.EchoTime ) && ( RepetitionTime == other.RepetitionTime ) && 
    ( AcquisitionNumber == other.AcquisitionNumber ) && ( this->GERawDataType == other.GERawDataType );
}

ImageFileDCM::ImageFileDCM( const char* filename )
{
  if ( cmtk::FileFormat::Identify( filename ) != cmtk::FILEFORMAT_DICOM )
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
    std::cerr << "Error: cannot read DICOM file " << filename << " (" << status.text() << ")" << std::endl;
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
  
  DcmStack stack;
  DcmTagKey searchKey;
    
  const char* tmpStr = NULL;
  if ( document->getValue( DCM_SeriesInstanceUID, tmpStr ) )
    SeriesUID = tmpStr;

  if ( document->getValue( DCM_SeriesDescription, tmpStr ) )
    SeriesDescription = tmpStr;

  if ( document->getValue( DCM_StudyInstanceUID, tmpStr ) )
    StudyUID = tmpStr;

  if ( document->getValue( DCM_EchoTime, tmpStr ) )
    EchoTime = tmpStr;

  if ( document->getValue( DCM_RepetitionTime, tmpStr ) )
    RepetitionTime = tmpStr;

  if ( document->getValue( DCM_ImagePositionPatient, tmpStr ) )
    ImagePositionPatient = tmpStr;

  if ( document->getValue( DCM_ImageOrientationPatient, tmpStr ) )
    ImageOrientationPatient = tmpStr;

  if ( ! document->getValue( DCM_InstanceNumber, ImageNumber ) )
    ImageNumber = 0;

  if ( ! document->getValue( DCM_AcquisitionNumber, AcquisitionNumber ) )
    AcquisitionNumber = 0;

  if ( WithExtensionsGE )
    {
    if ( ! document->getValue( DCM_RawDataType_ImageType, this->GERawDataType ) )
      this->GERawDataType = 3; // assume this is a magnitude image
    this->GERawDataType = std::min( 3, std::max( 0, (int)GERawDataType ) );
    }
  else
    {
    // GE extensions disabled
    this->GERawDataType = 0;
    }
}

ImageFileDCM::~ImageFileDCM()
{
  if ( fname ) free( fname );
  if ( fpath ) free( fpath );
}

class VolumeDCM : public std::vector<ImageFileDCM*> 
{
public:
  void AddImageFileDCM( ImageFileDCM *const image );
  
  bool Match ( const ImageFileDCM *newImage ) const;

  void WriteToArchive ( const std::string& name ) const;

  void print() const;
};

bool
VolumeDCM::Match ( const ImageFileDCM *newImage ) const
{
  if ( empty() ) 
    return 1;

  const ImageFileDCM *check = front();
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
VolumeDCM::AddImageFileDCM ( ImageFileDCM *const newImage )
{
  iterator it = begin();
  for ( ; it != end(); ++it )
    if ( newImage->ImageNumber < (*it)->ImageNumber ) break;
  insert( it, newImage );
}

void
VolumeDCM::WriteToArchive( const std::string& fname ) const
{
  cmtk::StudyImageSet studyImageSet;

  const ImageFileDCM *first = this->front();

  studyImageSet.SetImageFormat( cmtk::FILEFORMAT_DICOM );
  studyImageSet.SetImageDirectory( first->fpath );
  studyImageSet.SetMultiFile( true );

  for ( const_iterator it = this->begin(); it != this->end(); ++it ) 
    {
    studyImageSet.push_back( (*it)->fname );
    }

  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeFromStudy::Read( &studyImageSet ) );
  if ( volume )
    {
    cmtk::VolumeIO::Write( volume, fname.c_str() );
    if ( Verbose )
      {
      cmtk::StdOut.printf( "\nOutput file:%s\nImage size: %3dx%3dx%3d pixels\nPixel size: %.4fx%.4fx%.4f mm\n\n", 
			   fname.c_str(), volume->m_Dims[0], volume->m_Dims[1], volume->m_Dims[2], volume->m_Delta[0], volume->m_Delta[1], volume->m_Delta[2] );
      }
    }
  else
    {
    cmtk::StdOut << "No valid volume was read.\n";
    }
  
  if ( Verbose )
    {
    cmtk::StdOut << "DICOM Information: \n";
    cmtk::StdOut << "  Series:      " << first->SeriesUID << "\n";
    cmtk::StdOut << "  Study:       " << first->StudyUID << "\n";
    cmtk::StdOut << "  Acquisition: " << first->AcquisitionNumber << "\n";
    cmtk::StdOut << "  TR / TE:     " << first->RepetitionTime << "ms /" << first->EchoTime << "ms\n";
    cmtk::StdOut << "  Position:    " << first->ImagePositionPatient << "\n";
    cmtk::StdOut << "  Orientation: " << first->ImageOrientationPatient << "\n";

    if ( WithExtensionsGE )
      {
      cmtk::StdOut << "  GE Raw Data Type: " << first->GERawDataType << "\n";
      }
    
    cmtk::StdOut << "\nImage List:\n";
    for ( const_iterator it = this->begin(); it != this->end(); ++it ) 
      {
      cmtk::StdOut << (*it)->fname << " ";
      }
    cmtk::StdOut << "\n====================================================\n";
    }
}

class VolumeList : 
  public std::vector<VolumeDCM*> 
{
public:
  void AddImageFileDCM( ImageFileDCM *const newImage );
  
  void WriteToArchive();
};

inline std::string &
replacein(std::string &s, const std::string &sub, const std::string &other)
{
  assert(!sub.empty());
  size_t b = 0;
  for (;;)
    {
    b = s.find(sub, b);
    if (b == s.npos) break;
    s.replace(b, sub.size(), other);
    b += other.size();
    }
  return s;
}

void
VolumeList::WriteToArchive() 
{
  int idx = 1;
  std::map< std::string,std::vector<const VolumeDCM*> > pathToVolumeMap;
  for ( const_iterator it = begin(); it != end(); ++it ) 
    {
    if ( (*it)->size() > 1 )
      {
      // replace place holders
      std::string path( OutPathPattern );
      replacein( path, "%D", (*it)[0][0]->SeriesDescription );
      replacein( path, "%R", (*it)[0][0]->RepetitionTime );
      replacein( path, "%E", (*it)[0][0]->EchoTime );
      replacein( path, "%T", GERawDataTypeString[(*it)[0][0]->GERawDataType] );
      // finally, replace non-path characters
      replacein( path, " ", "_" );      
      replacein( path, ":", "_" );
      
      if ( path.length() > PATH_MAX )
	cmtk::StdErr << "ERROR: output path exceeds maximum path length";
      else
	pathToVolumeMap[path].push_back( *it );
      }
    }
  
  for ( std::map< std::string,std::vector<const VolumeDCM*> >::const_iterator it = pathToVolumeMap.begin(); it != pathToVolumeMap.end(); ++it )
    {						
    const size_t nVolumes = it->second.size();

    // if there is only one volume with the given output path, just write it
    if ( nVolumes == 1 )
      {					       
      // if there's a "number" tag, get rid of it.
      std::string uniquePath = it->first;
      replacein( uniquePath, "%N", "" );

      char finalPath[PATH_MAX];
      sprintf( finalPath, uniquePath.c_str(), idx++ );
      it->second[0]->WriteToArchive( finalPath );
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
	it->second[i]->WriteToArchive( finalPath );
	}
      }
    }
}


void
VolumeList::AddImageFileDCM( ImageFileDCM *const newImage )
{
  if ( empty() ) 
    {
    VolumeDCM *newVolumeDCM = new VolumeDCM;
    newVolumeDCM->AddImageFileDCM( newImage );
    push_back( newVolumeDCM );
    } 
  else
    {
    const_iterator it = begin();
    while ( it != end() ) 
      {
      VolumeDCM *study = *it;
      if ( study->Match( newImage ) ) 
	{
	study->AddImageFileDCM( newImage );
	return;
	} 
      else 
	{
	++it;
	}
      }
    VolumeDCM *newVolumeDCM = new VolumeDCM;
    newVolumeDCM->AddImageFileDCM( newImage );
    push_back( newVolumeDCM );    
    }
}

int
traverse_directory( VolumeList& studylist, const char *path, const char *wildcard )
{
  char fullname[PATH_MAX];

  std::vector<std::string> fileNameList;

#ifdef _MSC_VER
  WIN32_FIND_DATA fData;
  char pattern[PATH_MAX];
  snprintf( pattern, sizeof( pattern ), "%s\\%s", path, wildcard );
  HANDLE hFind = FindFirstFile( pattern, &fData);
  do
    {
    snprintf( fullname, sizeof( fullname ), "%s\\%s", path, fData.cFileName );
    if ( fData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY )
      {
      if ( Recursive && (fData.cFileName[0] != '.') )
	{
	traverse_directory( studylist, fullname, wildcard );
	}
      }
    else
      {
      fileNameList.push_back( fullname );
      (cmtk::StdErr << "\r" << progress_chars[ ++progress % 4 ]).flush();
      }
    }
  while (FindNextFile(hFind, &fData) != 0);
#else    
  DIR *dir_pointer = opendir ( path );
  if ( dir_pointer != NULL ) 
    {
    struct dirent *entry_pointer;

    while ( (entry_pointer = readdir(dir_pointer)) ) 
      {
      strcat( strcat( strcpy( fullname, path ), "/"), entry_pointer->d_name );
      struct stat entry_status;
      if ( !stat(fullname, &entry_status) ) 
	{
	if ( S_ISDIR( entry_status.st_mode ) && Recursive && (entry_pointer->d_name[0] != '.') ) 
	  {
	  strcat( fullname, "/" );
	  traverse_directory( studylist, fullname, wildcard );
	  } 
	else
	  {
	  if ( !fnmatch(wildcard,entry_pointer->d_name,FNM_PERIOD) ) 
	    {
	    cmtk::StdErr << "\r" << progress_chars[ ++progress % 4 ];
	    cmtk::StdErr.flush();
	    fileNameList.push_back( fullname );
	    }
	  }
	}
      }
    (void) closedir(dir_pointer);
    }
#endif

  if ( SortFilesByName )
    {
    std::sort( fileNameList.begin(), fileNameList.end() );
    }
  for ( std::vector<std::string>::const_iterator it = fileNameList.begin(); it != fileNameList.end(); ++it )
    {
    try 
      {
      studylist.AddImageFileDCM( new ImageFileDCM( it->c_str() ) );
      }
    catch (int)
      {
      }
    }
  std::cout << "\r";
  return 0;
}


int
main ( int argc, char *argv[] )
{
  if (! dcmDataDict.isDictionaryLoaded() ) 
    {
#ifdef CMAKE_DCMDICTPATH
    if ( dcmDataDict.wrlock().loadDictionary( CMAKE_DCMDICTPATH ) )
      {
      dcmDataDict.unlock();
      }
    else
#endif
      {
      std::cerr << "Data dictionary not avaliable. Please set DCMDICTPATH variable as path to dicom.dic file.\n";
      exit( 1 );
      }
    }

  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "DICOM to Image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Combine sets of DICOM slices to 3D image stacks" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] directory" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddOption( Key( 'O', "out-pattern" ), &OutPathPattern, "Output image path pattern (using printf substitutions) [default: '%03d.hdr']" );
    cl.AddSwitch( Key( 'r', "recurse" ), &Recursive, true, "Recurse into directories [default: no]" );
    cl.AddSwitch( Key( 'n', "no-sort" ), &SortFilesByName, false, "Do NOT sort files by file name (determines order when resolving spatial collisions)" );
    cl.AddSwitch( Key( 's', "sort" ), &SortFilesByName, true, "Sort files by file name (this is the default)" );
    
    cl.AddSwitch( Key( "no-orientation-check" ), &DisableOrientationCheck, true, "Disable checking of image orientations (to avoid rounding issues)" );
    cl.AddOption( Key( "tolerance" ), &Tolerance, "Tolerance for floating-point comparisons (must be >= 0; 0 = exact matches only; default: 1e-5)" );
    cl.AddSwitch( Key( "ge-extensions" ), &WithExtensionsGE, true, "Enable GE extensions (e.g., detect image type magnitude vs. complex)" );

    if ( ! cl.Parse() )
      return 1;

    const char* next = cl.GetNext();
    if ( next )
      {
      SearchRootDir = next;
      }
    }
  catch ( cmtk::CommandLine::Exception& ex )
    {
    cmtk::StdErr << ex << "\n";
    exit( 1 );
    }
  
  if ( WithExtensionsGE && !dcmDataDict.rdlock().findEntry( "RawDataType_ImageType" ) )
    {
    dcmDataDict.unlock();
    dcmDataDict.wrlock().addEntry( new DcmDictEntry( 0x0043, 0x102f, EVR_SS, "RawDataType_ImageType", 1, 1, NULL, OFFalse, "GE" ) );
    dcmDataDict.unlock();
    }

  VolumeList studylist;
  
  traverse_directory( studylist, SearchRootDir, "*" );
  
  studylist.WriteToArchive();
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace dcm2image
} // namespace apps
} // namespace cmtk
#endif
