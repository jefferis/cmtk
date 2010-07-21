/*
//
//  Copyright 2004-2010 SRI International
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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include "cmtkVolumeFromStudy.h"

#include "IO/cmtkMountPoints.h"
#include "IO/cmtkVolumeIO.h"
#include "IO/cmtkVolumeFromFile.h"
#include "IO/cmtkDICOM.h"

#include "System/cmtkCompressedStream.h"
#include "System/cmtkProgress.h"
#include "System/cmtkConsole.h"

#include "Base/cmtkVolume.h"
#include "Base/cmtkLandmarkList.h"
#include "Base/cmtkUniformVolume.h"

#include <cstring>
#include <climits>

#ifdef DEBUG
#  include <cstdio>
#endif

#include <memory>

namespace
cmtk
{

/** \addtogroup IO */
//@{

const UniformVolume::SmartPtr
VolumeFromStudy::Read
( const Study* study, const bool verbose )
{
  if ( !study ) 
    return UniformVolume::SmartPtr( NULL );

  const StudyImageSet* studyImageSet = dynamic_cast<const StudyImageSet*>( study );
  if ( studyImageSet ) 
    {
    VolumeFromStudy vfs;    
    return vfs.AssembleVolume( studyImageSet );
    } 
  else
    return VolumeIO::Read( study->GetFileSystemPath(), verbose );
}

const UniformVolume::SmartPtr
VolumeFromStudy::AssembleVolume( const StudyImageSet* study, const bool verbose )
{
  UniformVolume::SmartPtr Result( NULL );
  
  if ( study->size() < 2 ) 
    return Result;
  
  try
    {
    if ( verbose )
      fprintf( stderr, "\rReading images from path %s ...\n", MountPoints::Translate( study->GetImageDirectory() ) );
    
    Progress::Begin( 0, study->size(), 1, "Volume image assembly" );
    
    unsigned int nextPlane = 0;
    StudyImageSet::const_iterator it = study->begin();
    while ( it != study->end() ) 
      {      
      if ( verbose )
	fprintf( stderr, "\r%s", it->c_str() );
      
      char fullpath[PATH_MAX];
      snprintf( fullpath, sizeof( fullpath ), "%s/%s", MountPoints::Translate( study->GetImageDirectory() ), it->c_str() );
      
      ScalarImage::SmartPtr image = ScalarImage::SmartPtr( DICOM::Read( fullpath ) );

      // TODO: when returning NULL here, we also should tell
      // VolumeFromSlices that we give up, so it can free its
      // temporary storage.
      if ( !image ) 
	return UniformVolume::SmartPtr( NULL );

      if ( ! nextPlane ) 
	{
	// special treatment for first image in sequence.
	if ( study->GetMultiFile() )
	  InitSequence( image, study->size() );
	else
	  InitSequence( image, study->m_Dims[AXIS_Z] );
	}
      
      const char *error = FillPlane( nextPlane, image );
      
      Progress::SetProgress( nextPlane );
      
      if ( error ) 
	{
	StdErr.printf( "ERROR: %s: %s\n", fullpath, error );
	return UniformVolume::SmartPtr( NULL );
	}
      
      ++it;
      }
    Progress::Done();
    
    Result = this->FinishVolume();
    
    TypedArray::SmartPtr data = Result->GetData();
    if ( data ) 
      {
      if ( study->GetPadding() && ! data->GetPaddingFlag() ) 
	{
	data->SetPaddingValue( study->GetPaddingValue() );
	}
      }
    }
  catch (...) 
    {
    }

  return Result;
}

} // namespace cmtk
