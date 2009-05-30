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

#include <cmtkconfig.h>

#include <cmtkTypedArray.h>
#include <cmtkTemplateArray.h>
#include <cmtkVector3D.h>

#include <cmtkTypedStreamStudylist.h>
#include <cmtkVolumeIO.h>

#include <cmtkXform.h>
#include <cmtkAffineXform.h>
#include <cmtkWarpXform.h>

#include <cmtkReformatVolume.h>
#include <cmtkUniformVolumeInterpolator.h>
#include <cmtkSincInterpolator.h>
#include <cmtkLinearInterpolator.h>
#include <cmtkCubicInterpolator.h>

#include <stdio.h>
#include <string.h>

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#include <iostream>
#include <list>
#include <vector>

#include <cmtkCommandLine.h>
#include <cmtkConsole.h>
#include <cmtkStrUtility.h>
#include <cmtkMemory.h>

bool Verbose = false;

bool IgnoreReference = false;
bool Normalize = false;

const char* SetReference = NULL;
cmtk::ScalarDataType DataType = cmtk::TYPE_NONE;

cmtk::Interpolators::InterpolationEnum Interpolation = cmtk::Interpolators::LINEAR;
const char *outImage = "average_grey.hdr";

int FromSlice = -1;
int ToSlice = -1;

const char* InListGlobal = NULL;
std::list<const char*> InListNames;

const char* ReplaceFrom;
std::map<std::string,std::string> ReplaceMap;

const char* AddReplaceFrom( const char* arg )
{
  ReplaceFrom = arg;
  return NULL;
}

const char* AddReplaceTo( const char* arg )
{
  ReplaceMap[std::string(ReplaceFrom)] = std::string( arg );
  return NULL;
}

void
GetNormalizationCoefficients
( cmtk::TypedArray::SmartPtr& floatingData, 
  const cmtk::Types::DataItem refMean, const cmtk::Types::DataItem refVariance,
  cmtk::Types::DataItem& scale, cmtk::Types::DataItem& offset )
{
  cmtk::Types::DataItem fltMean, fltVariance;
  floatingData->GetStatistics( fltMean, fltVariance );
  scale = sqrt( refVariance ) / sqrt( fltVariance );
  offset = refMean - scale * fltMean;

  if ( Verbose ) 
    {
    fprintf( stderr, "Converted grey values: %f +- %f -> %f +- %f\n", fltMean, fltVariance, refMean, refVariance );
    }
}

void
ReformatAndAdd
( cmtk::UniformVolume::SmartPtr& referenceVolume, 
  cmtk::UniformVolume::SmartPtr& floatingVolume, 
  cmtk::Xform::SmartPtr& xform, float *const averageData, byte *const countData,
  const cmtk::Types::DataItem normFactor = 1.0, const cmtk::Types::DataItem normOffset = 0.0 )
{
  cmtk::Types::DataItem value;
  cmtk::TypedArray::SmartPtr floatingData = floatingVolume->GetData();

  const cmtk::UniformVolumeInterpolatorBase::SmartPtr interpolator( cmtk::ReformatVolume::CreateInterpolator( Interpolation, floatingVolume ) );

  if ( xform ) 
    {
    const int *dims = referenceVolume->GetDims();
    xform->RegisterVolume( referenceVolume );

#pragma omp parallel for private(value)
    for ( int pZ = 0; pZ < dims[cmtk::AXIS_Z]; ++pZ ) 
      {
      std::vector<cmtk::Vector3D> pFlt( dims[cmtk::AXIS_X] );
      size_t offset = pZ * dims[cmtk::AXIS_X] * dims[cmtk::AXIS_Y];
      for ( int pY = 0; pY < dims[cmtk::AXIS_Y]; ++pY ) 
	{
	xform->GetTransformedGridSequence( &(pFlt[0]), dims[cmtk::AXIS_X], 0, pY, pZ );
	
	for ( int pX = 0; pX < dims[cmtk::AXIS_X]; ++pX, ++offset ) 
	  {
	  if ( interpolator->GetDataAt( pFlt[pX], value ) )
	    {
	    averageData[offset] += ((value * normFactor) + normOffset);
	    ++countData[offset];
	    }
	  }
	}
      }
    }
  else
    {
#pragma omp parallel for private(value)
    for ( size_t index = 0; index < referenceVolume->GetNumberOfPixels(); ++ index ) 
      {
      if ( floatingData->Get( value, index ) ) 
	{
	averageData[index] += value;
	++countData[index];
	}
      }
    } 
}    

void
ReformatAndAdd
( cmtk::UniformVolume::SmartPtr& referenceVolume, 
  cmtk::Xform::SmartPtr& xformGlobal,
  cmtk::UniformVolume::SmartPtr& floatingVolume, 
  cmtk::Xform::SmartPtr& xform, float *const averageData, byte *const countData,
  const cmtk::Types::DataItem normFactor = 1.0, const cmtk::Types::DataItem normOffset = 0.0 )
{
  cmtk::Types::DataItem value;
  cmtk::TypedArray::SmartPtr floatingData = floatingVolume->GetData();

  const cmtk::UniformVolumeInterpolatorBase::SmartPtr interpolator( cmtk::ReformatVolume::CreateInterpolator( Interpolation, floatingVolume ) );  
  if ( xform ) 
    {
    const int *dims = referenceVolume->GetDims();

#pragma omp parallel for private(value)
    for ( int pZ = 0; pZ < dims[cmtk::AXIS_Z]; ++pZ ) 
      {
      cmtk::Vector3D pFlt;
      size_t offset = pZ * dims[cmtk::AXIS_X] * dims[cmtk::AXIS_Y];
      for ( int pY = 0; pY < dims[cmtk::AXIS_Y]; ++pY ) 
	for ( int pX = 0; pX < dims[cmtk::AXIS_X]; ++pX, ++offset ) 
	  {
	  referenceVolume->GetGridLocation( pFlt, pX, pY, pZ );
	  xformGlobal->ApplyInPlace( pFlt );
	  
	  bool success = xform->InDomain( pFlt );
	  if ( success )
	    {
	    xform->ApplyInPlace( pFlt );
	    success = interpolator->GetDataAt( pFlt, value );
	    }
	  if ( success ) 
	    {
	    averageData[offset] += ((value * normFactor) + normOffset);
	    ++countData[offset];
	    }
	  }
      }
    }
}    

int
main ( const int argc, const char* argv[] ) 
{
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Average grey-level images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] list0 [list1 ...]" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode." );

    cl.AddOption( Key( 'o', "output" ), &outImage, "Output image file name" );

    cl.AddOption( Key( 'f', "from-slice" ), &FromSlice, "Work _from_ slice index" );
    cl.AddOption( Key( 't', "to-slice" ), &FromSlice, "Work _to_ slice index" );
    cl.AddOption( Key( 'r', "set-reference" ), &SetReference, "Set reference grid to file <image>." );
    cl.AddOption( Key( 'X', "global-xform" ), &InListGlobal, "Set global transformation applied to all images prior to individual transformation." );

    cl.AddSwitch( Key( 'c', "cubic" ), &Interpolation, cmtk::Interpolators::CUBIC, "Use cubic interpolation" );
    cl.AddSwitch( Key( 'n', "normalize" ), &Normalize, true, "Normalize intensity ranges accross images." );
    cl.AddSwitch( Key( 'i', "ignore-reference" ), &IgnoreReference, true, "Ignore reference image data. Only average floating images." );

    cl.AddSwitch( Key( 'C', "char" ), &DataType, cmtk::TYPE_CHAR, "8 bits, signed" );
    cl.AddSwitch( Key( 'B', "byte" ), &DataType, cmtk::TYPE_BYTE, "8 bits, unsigned" );
    cl.AddSwitch( Key( 'S', "short" ), &DataType, cmtk::TYPE_SHORT, "16 bits, signed" );
    cl.AddSwitch( Key( 'U', "ushort" ), &DataType, cmtk::TYPE_USHORT, "16 bits, unsigned" );
    cl.AddSwitch( Key( 'I', "int" ), &DataType, cmtk::TYPE_INT, "32 bits signed" );
    cl.AddSwitch( Key( 'F', "float" ), &DataType, cmtk::TYPE_FLOAT, "32 bits floating point" );
    cl.AddSwitch( Key( 'D', "double" ), &DataType, cmtk::TYPE_DOUBLE, "64 bits floating point\n" );

    cl.AddCallback( Key( "replace-from" ), AddReplaceFrom, "Replace from pattern" );
    cl.AddCallback( Key( "replace-to" ), AddReplaceTo, "Replace to pattern" );

    cl.Parse();

    const char *next = cl.GetNext();
    while ( next )
      {
      InListNames.push_back( next );
      next = cl.GetNextOptional();
      }
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e;
    exit( 1 );
    }

  cmtk::TypedStreamStudylist **studylist = cmtk::Memory::AllocateArray<cmtk::TypedStreamStudylist*>( InListNames.size() );
  
  const char* referenceStudy = NULL;
  size_t idx = 0;
  for ( std::list<const char*>::const_iterator it = InListNames.begin(); it != InListNames.end(); ++it, ++idx )
    {
    if ( Verbose ) 
      std::cerr << "Opening studylist " << *it << ".\n";
    
    studylist[idx] = new cmtk::TypedStreamStudylist();
    studylist[idx]->Read( *it );

    if ( ! idx ) 
      {
      referenceStudy = studylist[0]->GetReferenceStudyPath();
      } 
    else
      {
      if ( strcmp( studylist[idx]->GetReferenceStudyPath(), referenceStudy ) ) 
	{
	std::cerr << "Studylist #" << idx << " has a different reference study.\n";
	exit( 1 );
	}
      }
    }
  
  cmtk::Xform::SmartPtr xformGlobal( NULL );
  if ( InListGlobal )
    {
    cmtk::TypedStreamStudylist globalStudyList( InListGlobal );
    xformGlobal = globalStudyList.GetWarpXform();
    if ( ! xformGlobal )
      xformGlobal = globalStudyList.GetAffineXform();
    }

  std::string actualPath;
  if ( SetReference )
    actualPath = SetReference;
  else
    actualPath = cmtk::StrReplace( referenceStudy, ReplaceMap );
  cmtk::UniformVolume::SmartPtr referenceVolume( cmtk::VolumeIO::ReadOriented( actualPath.c_str(), Verbose ) );
  if ( ! referenceVolume ) 
    {
    std::cerr << "Could not read reference image " << referenceStudy << "\n";
    exit( 1 );
    }

  size_t reformattedPixelCount = referenceVolume->GetNumberOfPixels();
  cmtk::SmartPointer<cmtk::FloatArray> averageArray
    ( dynamic_cast<cmtk::FloatArray*>( cmtk::TypedArray::Create( cmtk::TYPE_FLOAT, reformattedPixelCount ) ) );
  averageArray->ClearArray( 0 );
  float* averageData = averageArray->GetDataPtrTemplate();
  
  byte* countData = cmtk::Memory::AllocateArray<byte>( reformattedPixelCount );
  memset( countData, 0, reformattedPixelCount );
  
  cmtk::TypedArray::SmartPtr referenceData = referenceVolume->GetData();
  cmtk::Types::DataItem refMean, refVariance;
  if ( Normalize )
    {
    if ( referenceData )
      {
      referenceData->GetStatistics( refMean, refVariance );
      }
    }
  
  if ( (DataType == cmtk::TYPE_NONE) && referenceData )
    {
    DataType = referenceData->GetType();
    }

  if ( referenceData && ! IgnoreReference ) 
    {
    cmtk::Xform::SmartPtr xform( NULL );
    ReformatAndAdd( referenceVolume, referenceVolume, xform, averageData, countData );
    }
  
  cmtk::Types::DataItem normFactor = 1.0, normOffset = 0;
  
  for ( size_t idx = 0; idx < InListNames.size(); ++idx ) 
    {
    std::string floatingName = cmtk::StrReplace( studylist[idx]->GetFloatingStudyPath(), ReplaceMap );
    cmtk::UniformVolume::SmartPtr floatingVolume( cmtk::VolumeIO::ReadOriented( floatingName.c_str(), Verbose ) );
    if ( floatingVolume && floatingVolume->GetData() ) 
      {
      if ( DataType == cmtk::TYPE_NONE )
	DataType = floatingVolume->GetData()->GetType();
      
      cmtk::Xform::SmartPtr xform = studylist[idx]->GetWarpXform();
      if ( ! xform )
	xform = cmtk::Xform::SmartPtr::DynamicCastFrom( studylist[idx]->GetAffineXform() );
      
      if ( Verbose ) 
	{
	std::cerr << "Reformating floating data #" << idx << "\r";
	}
      
      if ( Normalize ) 
	{
	cmtk::TypedArray::SmartPtr floatingData = floatingVolume->GetData();
	if ( ! referenceData && ! idx )
	  {
	  floatingData->GetStatistics( refMean, refVariance );
	  }
	GetNormalizationCoefficients( floatingData, refMean, refVariance, normFactor, normOffset );
	}

      if ( xformGlobal )
	ReformatAndAdd( referenceVolume, xformGlobal, floatingVolume, xform, averageData, countData, normFactor, normOffset );
      else
	ReformatAndAdd( referenceVolume, floatingVolume, xform, averageData, countData, normFactor, normOffset );
      } 
    else
      {
      std::cerr << "ERROR: Could not read floating volume " << floatingName << "\n";
      exit( 1 );
      }
    }
  
  float *avgPtr = averageData;
  byte *countPtr = countData;
  for ( unsigned int idx = 0; idx < reformattedPixelCount; ++idx, ++avgPtr, ++countPtr ) 
    {
    if ( *countPtr ) 
      {
      (*avgPtr) /= (1.0 * (*countPtr) );
      }
    }
  
  cmtk::Types::DataItem min, max;
  averageArray->GetRange( min, max );
  //  cmtk::Types::DataItem factor = 255.0 / (max-min);
  //  averageArray->Rescale( factor, -min / factor, 0, 255 );

  cmtk::TypedArray::SmartPtr outputData( averageArray->Convert( DataType ) );
  referenceVolume->SetData( outputData );
  
  cmtk::VolumeIO::Write( referenceVolume, outImage, Verbose );

  for ( size_t idx = 0; idx < InListNames.size(); ++idx )
    delete studylist[idx];
  delete[] studylist;
  
  delete[] countData;
}

