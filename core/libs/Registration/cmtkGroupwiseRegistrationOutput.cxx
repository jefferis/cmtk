/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include "cmtkGroupwiseRegistrationOutput.h"

#include <System/cmtkDebugOutput.h>

#include <IO/cmtkGroupwiseRegistrationFunctionalIO.h>
#include <IO/cmtkClassStream.h>
#include <IO/cmtkStudyList.h>
#include <IO/cmtkClassStreamStudyList.h>
#include <IO/cmtkVolumeIO.h>

#include <limits.h>
#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

bool
GroupwiseRegistrationOutput::WriteGroupwiseArchive( const char* path ) const
{
#ifdef CMTK_USE_MPI    
  // only root process needs to write outputs
  if ( MPI::COMM_WORLD.Get_rank() == 0 )
#endif
    {
    // create class stream archive.
    if ( path )
      {
      ClassStream stream;

      if ( this->m_OutputRootDirectory )
	{
	char completePath[PATH_MAX];
	snprintf( completePath, sizeof( completePath ), "%s%c%s", this->m_OutputRootDirectory, (int)CMTK_PATH_SEPARATOR, path );
	stream.Open( completePath, ClassStream::WRITE );
	}
      else
	stream.Open( path, ClassStream::WRITE );
      
      if ( ! stream.IsValid() ) return false;
      stream << *this->m_Functional;
      stream.Close();
      }
    }
  return true;
}

bool 
GroupwiseRegistrationOutput::WriteXformsSeparateArchives
( const char* path, const char* templatePath )
{ 
#ifdef CMTK_USE_MPI    
  // only root process needs to write outputs
  if ( MPI::COMM_WORLD.Get_rank() == 0 )
#endif
    {
    if ( path )
      {
      char fullPath[PATH_MAX];
      
      for ( size_t img = 0; img < this->m_Functional->GetNumberOfTargetImages(); ++img )
	{
	StudyList slist;
	Study::SmartPtr refstudy;
	if ( this->m_OutputRootDirectory  && ! this->m_ExistingTemplatePath )
	  {
	  snprintf( fullPath, sizeof( fullPath ), "%s%c%s", this->m_OutputRootDirectory, CMTK_PATH_SEPARATOR, templatePath );
	  refstudy = slist.AddStudy( fullPath );
	  }
	else
	  {
	  refstudy = slist.AddStudy( templatePath );
	  }
	
	const UniformVolume* image = this->m_Functional->GetOriginalTargetImage( img );
	Study::SmartPtr imgstudy = slist.AddStudy( image->GetMetaInfo( META_FS_PATH ).c_str() );
	
	WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom( this->m_Functional->GetGenericXformByIndex( img ) );
	if ( warpXform )
	  {
	  AffineXform::SmartPtr affineXform( warpXform->GetInitialAffineXform() );
	  slist.AddXform( refstudy, imgstudy, affineXform, warpXform );
	  }
	else
	  {
	  AffineXform::SmartPtr affineXform = AffineXform::SmartPtr::DynamicCastFrom( this->m_Functional->GetGenericXformByIndex( img ) );
	  slist.AddXform( refstudy, imgstudy, affineXform );
	  }
	
	if ( this->m_OutputRootDirectory )
	  {
	  snprintf( fullPath, sizeof( fullPath ), "%s%c%s%ctarget-%03d.list", this->m_OutputRootDirectory, CMTK_PATH_SEPARATOR, path, CMTK_PATH_SEPARATOR, (int)img );
	  }
	else
	  {
	  snprintf( fullPath, sizeof( fullPath ), "%s%ctarget-%03d.list", path, CMTK_PATH_SEPARATOR, (int)img );
	  }
	ClassStreamStudyList::Write( fullPath, &slist );
	}
      }
    }
  return true;
}

bool
GroupwiseRegistrationOutput::WriteAverageImage( const char* path, const cmtk::Interpolators::InterpolationEnum interp, const bool useTemplateData )
{
  // reformat output and generate average images
  if ( path )
    {
    UniformVolume::SmartPtr templateGrid = this->m_Functional->GetTemplateGrid();
    const size_t numberOfPixels = templateGrid->GetNumberOfPixels();

    TypedArray::SmartPtr average( TypedArray::Create( TYPE_FLOAT, numberOfPixels ) );
    float* averagePtr = static_cast<float*>( average->GetDataPtr() );

    TypedArray::SmartPtr count( TypedArray::Create( TYPE_USHORT, numberOfPixels ) );
    unsigned short* countPtr = static_cast<unsigned short*>( count->GetDataPtr() );
    
    if ( useTemplateData )
      {
      if ( ! templateGrid->GetData() )
	{
	UniformVolume::SmartPtr readImage( VolumeIO::ReadOriented( templateGrid->GetMetaInfo( META_FS_PATH ).c_str() ) );
	templateGrid->SetData( readImage->GetData() );
	}

      for ( size_t px = 0; px < numberOfPixels; ++px )
	{
	averagePtr[px] = static_cast<float>( templateGrid->GetDataAt( px ) );
	}
      count->Fill( 1 );
      }  
    else
      {
      average->Fill( 0 );
      count->Fill( 0 );
      }

    DebugOutput( 1 ) << "Reformating output images\n";

#ifdef CMTK_USE_MPI
    const size_t idxFrom = MPI::COMM_WORLD.Get_rank();
    const size_t idxSkip = MPI::COMM_WORLD.Get_size();
#else
    const size_t idxFrom = 0;
    const size_t idxSkip = 1;
#endif
    for ( size_t idx = idxFrom; idx < this->m_Functional->GetNumberOfTargetImages(); idx += idxSkip )
      {
      UniformVolume::SmartPtr floatingVolume = this->m_Functional->GetOriginalTargetImage( idx );
      if ( !floatingVolume->GetData() )
	floatingVolume = UniformVolume::SmartPtr( VolumeIO::ReadOriented( floatingVolume->GetMetaInfo( META_FS_PATH ).c_str() ) );
      
      cmtk::ReformatVolume reformat;
      reformat.SetReferenceVolume( templateGrid );
      reformat.SetFloatingVolume( floatingVolume );
      reformat.SetInterpolation( interp );
      
      AffineXform::SmartPtr affineXform = AffineXform::SmartPtr::DynamicCastFrom( this->m_Functional->GetGenericXformByIndex( idx ) );
      if ( affineXform )
	{
	reformat.SetAffineXform( affineXform );
	}

      WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom( this->m_Functional->GetGenericXformByIndex( idx ) );
      if ( warpXform )
	reformat.SetWarpXform( warpXform );
      
      UniformVolume::SmartPtr ref( reformat.PlainReformat() );
      const TypedArray* data = ref->GetData();
#pragma omp parallel for
      for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
	{
	Types::DataItem v;
	if ( data->Get( v, i ) )
	  {
	  averagePtr[i] += static_cast<float>( v );
	  ++countPtr[i];
	  }
	}
      }

#ifdef CMTK_USE_MPI
    float* averagePtrMPI = Memory::ArrayC::Allocate<float>( numberOfPixels );
    MPI::COMM_WORLD.Reduce( averagePtr, averagePtrMPI, numberOfPixels, MPI::FLOAT, MPI::SUM, 0 );
    memcpy( averagePtr, averagePtrMPI, numberOfPixels * sizeof( *averagePtr ) );
    Memory::ArrayC::Delete( averagePtrMPI );

    unsigned short* countPtrMPI = Memory::ArrayC::Allocate<unsigned short>( numberOfPixels );
    MPI::COMM_WORLD.Reduce( countPtr, countPtrMPI, numberOfPixels, MPI::UNSIGNED_SHORT, MPI::SUM, 0 );
    memcpy( countPtr, countPtrMPI, numberOfPixels * sizeof( *countPtr ) );
    Memory::ArrayC::Delete( countPtrMPI );
#endif

#ifdef CMTK_USE_MPI
    if ( MPI::COMM_WORLD.Get_rank() == 0 )
#endif
      {
#pragma omp parallel for
      for ( int i = 0; i < static_cast<int>( numberOfPixels ); ++i )
	{
	if ( countPtr[i] )
	  averagePtr[i] /= countPtr[i];
	else
	  average->SetPaddingAt( i );
	}
      templateGrid->SetData( average );

      if ( this->m_OutputRootDirectory )
	{
	char fullPath[PATH_MAX];
	snprintf( fullPath, sizeof( fullPath ), "%s%c%s", this->m_OutputRootDirectory, CMTK_PATH_SEPARATOR, path );
	VolumeIO::Write( *templateGrid, fullPath );
	}
      else
	{
	VolumeIO::Write( *templateGrid, path );
	}
      }
    }
  
  return 0;
}

} // namespace cmtk
