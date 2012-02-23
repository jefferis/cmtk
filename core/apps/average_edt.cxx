/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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
#include <System/cmtkStrUtility.h>
#include <System/cmtkTimers.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkUniformDistanceMap.h>
#include <Base/cmtkTypedArray.h>
#include <Base/cmtkTemplateArray.h>
#include <Base/cmtkLinearInterpolator.h>
#include <Base/cmtkUniformVolumeInterpolator.h>
#include <Base/cmtkXformUniformVolume.h>
#include <Base/cmtkAffineXformUniformVolume.h>
#include <Base/cmtkSplineWarpXformUniformVolume.h>

#include <IO/cmtkFileFormat.h>
#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkTypedStreamStudylist.h>

#include <Registration/cmtkReformatVolume.h>

#include <math.h>
#include <list>

#ifdef CMTK_USE_GCD
#  include <dispatch/dispatch.h>
#endif

const char* DownsampleVolumeStr = NULL;
int DownsampleVolume[3] = { 1, 1, 1 };

bool InterpolateDistance = true;
float FeatureWindowRadius = -1;
bool ScaleIntensitiesToLabels = false;

const char *OutputFileName = NULL;

const char* WriteDistanceMapNameMask = NULL;
const char* WriteLabelMapNameMask = NULL;

std::list<const char*> InputFileList;

short NumberOfLabels = 0;

const char* ReplaceFrom;
std::map<std::string,std::string> ReplaceMap;

bool PaddingFlag = false;
float PaddingValue = 0;

void AddReplaceFrom( const char* arg )
{
  ReplaceFrom = arg;
}

void AddReplaceTo( const char* arg )
{
  ReplaceMap[std::string(ReplaceFrom)] = std::string( arg );
}

cmtk::TypedArray::SmartPtr
Average
( std::list<cmtk::UniformVolume::SmartPtr> volumes, const short numLabels )
{
  const int distanceMapFlags = cmtk::UniformDistanceMap<float>::VALUE_EXACT + cmtk::UniformDistanceMap<float>::SIGNED;

  const size_t numPixels = (*volumes.begin())->GetNumberOfPixels();

  bool labelFlags[256];
  memset( labelFlags, 0, sizeof( labelFlags ) );
  bool* pFlags = labelFlags; // need to use ptr to array to work around GCD bug

  for ( std::list<cmtk::UniformVolume::SmartPtr>::const_iterator it = volumes.begin(); it != volumes.end(); ++it )
    {
    const cmtk::TypedArray* data = (*it)->GetData();
#ifdef CMTK_USE_GCD
    const cmtk::Threads::Stride stride( numPixels );
    dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		    { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			for ( int i = 0; i < static_cast<int>( numPixels ); ++i )
#endif
      {
      cmtk::Types::DataItem l;
      if ( data->Get( l, i ) )
	pFlags[static_cast<short>( l )] = true;
      }
#ifdef CMTK_USE_GCD
		    });
#endif
    }

  cmtk::TypedArray::SmartPtr result( cmtk::TypedArray::Create( cmtk::TYPE_SHORT, numPixels ) );
  result->BlockSet( 0 /*value*/, 0 /*idx*/, numPixels /*len*/ );
  short* resultPtr = static_cast<short*>( result->GetDataPtr() );
  
  cmtk::FloatArray::SmartPtr totalDistance( new cmtk::FloatArray( numPixels ) );
  float* totalDistancePtr = totalDistance->GetDataPtrTemplate();

  cmtk::FloatArray::SmartPtr inOutDistance( new cmtk::FloatArray(numPixels ) );
  float* inOutDistancePtr = inOutDistance->GetDataPtrTemplate();

  totalDistance->BlockSet( 0 /*value*/, 0 /*idx*/, numPixels /*len*/ );
  for ( int label = 0; label < numLabels; ++label )
    {
    /// skip labels that are not in any image.
    if ( ! labelFlags[label] ) continue;

    cmtk::DebugOutput( 1 ) << "Processing label #" << label << "\r";

    inOutDistance->BlockSet( 0 /*value*/, 0 /*idx*/, numPixels /*len*/ );

    for ( std::list<cmtk::UniformVolume::SmartPtr>::const_iterator it = volumes.begin(); it != volumes.end(); ++it )
      {
      cmtk::UniformVolume::SmartPtr signedDistanceMap = cmtk::UniformDistanceMap<float>( *(*it), distanceMapFlags, label ).Get();
      const float* signedDistancePtr = (const float*)signedDistanceMap->GetData()->GetDataPtr();

      // if this is the first label, write directly to accumulation distance map
      if ( !label )
	{
#ifdef CMTK_USE_GCD
    const cmtk::Threads::Stride stride( numPixels );
    dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		    { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			for ( int i = 0; i < static_cast<int>( numPixels ); ++i )
#endif
	  {
	  totalDistancePtr[i] += signedDistancePtr[i];
	  }
#ifdef CMTK_USE_GCD
		    });
#endif
	}
      else
	// for all other labels, add to label distance map
	{
#ifdef CMTK_USE_GCD
    const cmtk::Threads::Stride stride( numPixels );
    dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		    { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			for ( int i = 0; i < static_cast<int>( numPixels ); ++i )
#endif
	  {
	  inOutDistancePtr[i] += signedDistancePtr[i];
	  }
#ifdef CMTK_USE_GCD
		    });
#endif
	}
      }

    // if this is not the first label, compare this label's sum distance map
    // (over all volumes) pixel by pixel and set this label where it is
    // closer than previous closest label
    if ( label )
      {
#ifdef CMTK_USE_GCD
      const cmtk::Threads::Stride stride( numPixels );
      dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		      { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			  for ( int i = 0; i < static_cast<int>( numPixels ); ++i )
#endif
	{
	if ( inOutDistancePtr[i] < totalDistancePtr[i] )
	  {
	  totalDistancePtr[i] = inOutDistancePtr[i];
	  resultPtr[i] = label;
	  }
	else
	  if ( !(inOutDistancePtr[i] > totalDistancePtr[i]) )
	    {
	    resultPtr[i] = numLabels;
	    }	  
	}
#ifdef CMTK_USE_GCD
		      });
#endif
      }

    if ( WriteDistanceMapNameMask )
      {
      cmtk::UniformVolume::SmartPtr distanceMap( (*volumes.begin())->CloneGrid() );
      distanceMap->SetData( totalDistance );

      char fname[PATH_MAX];
      if ( snprintf( fname, PATH_MAX, WriteDistanceMapNameMask, label ) > PATH_MAX )
	{
	cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
	}
      else
	{
	cmtk::VolumeIO::Write( *distanceMap, fname );
	}
      }

    if ( WriteLabelMapNameMask )
      {
      cmtk::UniformVolume::SmartPtr labelMap( (*volumes.begin())->CloneGrid() );
      labelMap->SetData( result );
      
      char fname[PATH_MAX];
      if ( snprintf( fname, PATH_MAX, WriteLabelMapNameMask, label ) > PATH_MAX )
	{
	cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
	}
      else
	{
	cmtk::VolumeIO::Write( *labelMap, fname );
	}
      }
    }
  
  return result;
}

cmtk::TypedArray::SmartPtr
AverageWindowed
( std::list<cmtk::UniformVolume::SmartPtr> volumes, const short numLabels )
{
  const int distanceMapFlags = cmtk::UniformDistanceMap<float>::VALUE_WINDOW + cmtk::UniformDistanceMap<float>::SIGNED;

  const size_t numPixels = (*volumes.begin())->GetNumberOfPixels();

  bool labelFlags[256];
  memset( labelFlags, 0, sizeof( labelFlags ) );
  bool *pFlags = labelFlags; // need pointer to work around GCD bug

  for ( std::list<cmtk::UniformVolume::SmartPtr>::const_iterator it = volumes.begin(); it != volumes.end(); ++it )
    {
    const cmtk::TypedArray* data = (*it)->GetData();

#ifdef CMTK_USE_GCD
    const cmtk::Threads::Stride stride( numPixels );
    dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		    { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			for ( int i = 0; i < static_cast<int>( numPixels ); ++i )
#endif
      {
      cmtk::Types::DataItem l;
      if ( data->Get( l, i ) )
	pFlags[static_cast<short>( l )] = true;
      }
#ifdef CMTK_USE_GCD
		    });
#endif
    }
  
  cmtk::TypedArray::SmartPtr result( cmtk::TypedArray::Create( cmtk::TYPE_FLOAT, numPixels ) );
  result->BlockSet( 0 /*value*/, 0 /*idx*/, numPixels /*len*/ );
  float* resultPtr = (float*) result->GetDataPtr();
  
  cmtk::TypedArray::SmartPtr resultDivider( cmtk::TypedArray::Create( cmtk::TYPE_SHORT, numPixels ) );
  resultDivider->BlockSet( 0 /*value*/, 0 /*idx*/, numPixels /*len*/ );
  short* resultDividerPtr = static_cast<short*>( resultDivider->GetDataPtr() );
  
  cmtk::FloatArray::SmartPtr totalDistance( new cmtk::FloatArray( numPixels ) );
  float* totalDistancePtr = totalDistance->GetDataPtrTemplate();
  
  cmtk::FloatArray::SmartPtr inOutDistance( new cmtk::FloatArray(numPixels ) );
  float* inOutDistancePtr = inOutDistance->GetDataPtrTemplate();
  
  totalDistance->BlockSet( 0 /*value*/, 0 /*idx*/, numPixels /*len*/ );
  for ( int label = 0; label < numLabels; ++label )
    {
    // skip labels that are not in any image.
    if ( ! labelFlags[label] ) continue;

    cmtk::DebugOutput( 1 ) << "Processing label #" << label << "\r";

    inOutDistance->BlockSet( 0 /*value*/, 0 /*idx*/, numPixels /*len*/ );

    std::list<cmtk::UniformVolume::SmartPtr>::const_iterator it;
    for ( it = volumes.begin(); it != volumes.end(); ++it )
      {
      cmtk::UniformVolume::SmartPtr signedDistanceMap = cmtk::UniformDistanceMap<float>( *(*it), distanceMapFlags, label, FeatureWindowRadius).Get();
      const float* signedDistancePtr = (const float*)signedDistanceMap->GetData()->GetDataPtr();

#ifdef CMTK_USE_GCD
      const cmtk::Threads::Stride stride( numPixels );
      dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		      { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			  for ( int i = 0; i < static_cast<int>( numPixels ); ++i )
#endif
			    {
			    inOutDistancePtr[i] += signedDistancePtr[i];
			    }
#ifdef CMTK_USE_GCD
		      });
#endif
      }
    
    // if this is not the first label, compare this label's sum distance map
    // (over all volumes) pixel by pixel and set this label where it is
    // closer than previous closest label
    if ( label )
      {
#ifdef CMTK_USE_GCD
      const cmtk::Threads::Stride stride( numPixels );
      dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		      { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			  for ( int i = 0; i < static_cast<int>( numPixels ); ++i )
#endif
	{	
	if ( inOutDistancePtr[i] < totalDistancePtr[i] )
	  {
	  totalDistancePtr[i] = inOutDistancePtr[i];
	  resultPtr[i] = label;
	  resultDividerPtr[i] = 1;
	  }
	else
	  {
	  if ( !(inOutDistancePtr[i] > totalDistancePtr[i]) )
	    {
	    resultPtr[i] += label;
	    ++resultDividerPtr[i];
	    }
	  }
	}
#ifdef CMTK_USE_GCD
		      });
#endif
      }

    if ( WriteDistanceMapNameMask )
      {
      cmtk::UniformVolume::SmartPtr distanceMap( (*volumes.begin())->CloneGrid() );
      distanceMap->SetData( totalDistance );

      char fname[PATH_MAX];
      if ( snprintf( fname, PATH_MAX, WriteDistanceMapNameMask, label ) > PATH_MAX )
	{
	cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
	}
      else
	{
	cmtk::VolumeIO::Write( *distanceMap, fname );
	}
      }

    if ( WriteLabelMapNameMask )
      {
      cmtk::UniformVolume::SmartPtr labelMap( (*volumes.begin())->CloneGrid() );
      labelMap->SetData( result );
      
      char fname[PATH_MAX];
      if ( snprintf( fname, PATH_MAX, WriteLabelMapNameMask, label ) > PATH_MAX )
	{
	cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
	}
      else
	{
	cmtk::VolumeIO::Write( *labelMap, fname );
	}
      }
    }

  // compute average between min and max window result
#ifdef CMTK_USE_GCD
  const cmtk::Threads::Stride stride( numPixels );
  dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		  { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
		      for ( int i = 0; i < static_cast<int>( numPixels ); ++i )
#endif
    {
    if ( resultDividerPtr[i] )
      {
      resultPtr[i] /= resultDividerPtr[i];
      }
    else
      {
      resultPtr[i] = numLabels;
      }
    resultDividerPtr[i] = 0;
    }
#ifdef CMTK_USE_GCD
		  });
#endif

  result = cmtk::TypedArray::SmartPtr( result->Convert( cmtk::TYPE_SHORT ) );
  return result;
}

cmtk::TypedArray::SmartPtr
Average
( const cmtk::UniformVolume* referenceVolume,
  std::list<cmtk::UniformVolume::SmartPtr> volumes,
  std::list<cmtk::XformUniformVolume::SmartConstPtr> xforms,
  const short numLabels )
{
  const int distanceMapFlags = cmtk::UniformDistanceMap<float>::VALUE_EXACT + cmtk::UniformDistanceMap<float>::SIGNED;

  const size_t nPixelsReference = referenceVolume->GetNumberOfPixels();
  const cmtk::DataGrid::IndexType& referenceDims = referenceVolume->GetDims();

  bool labelFlags[256];
  memset( labelFlags, 0, sizeof( labelFlags ) );
  bool *pFlags = labelFlags; // need pointer to work around GCD bug

  for ( std::list<cmtk::UniformVolume::SmartPtr>::const_iterator itV = volumes.begin(); itV != volumes.end(); ++itV )
    {
    const cmtk::TypedArray* data = (*itV)->GetData();

#ifdef CMTK_USE_GCD
    const cmtk::Threads::Stride stride( data->GetDataSize() );
    dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		    { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			for ( int i = 0; i < static_cast<int>( data->GetDataSize() ); ++i )
#endif
      {
      cmtk::Types::DataItem l;
      if ( data->Get( l, i ) )
	pFlags[static_cast<short>( l )] = true;
      }
#ifdef CMTK_USE_GCD
		    });
#endif
    }
  
  cmtk::TypedArray::SmartPtr result( cmtk::TypedArray::Create( cmtk::TYPE_SHORT, nPixelsReference ) );
  result->BlockSet( 0 /*value*/, 0 /*idx*/, nPixelsReference /*len*/ );
  short* resultPtr = static_cast<short*>( result->GetDataPtr() );
  
  cmtk::FloatArray::SmartPtr referenceInOutDistance( new cmtk::FloatArray( nPixelsReference ) );
  float* referenceInOutDistancePtr = referenceInOutDistance->GetDataPtrTemplate();

  cmtk::FloatArray::SmartPtr totalDistance ( new cmtk::FloatArray( nPixelsReference ) );
  float* totalDistancePtr = totalDistance->GetDataPtrTemplate();

  totalDistance->BlockSet( 0 /*value*/, 0 /*idx*/, nPixelsReference /*len*/ );
  for ( int label = 0; label < numLabels; ++label )
    {
    /// skip labels that are not in any image.
    if ( ! labelFlags[label] ) continue;

    cmtk::DebugOutput( 1 ) << "Processing label #" << label << "\r";

    referenceInOutDistance->BlockSet( 0 /*value*/, 0 /*idx*/, nPixelsReference /*len*/ );

    std::list<cmtk::XformUniformVolume::SmartConstPtr>::const_iterator itX = xforms.begin();
    for ( std::list<cmtk::UniformVolume::SmartPtr>::const_iterator itV = volumes.begin(); itV != volumes.end(); ++itV, ++itX )
      {
      cmtk::UniformVolume::SmartPtr signedDistanceMap = cmtk::UniformDistanceMap<float>( *(*itV), distanceMapFlags, label ).Get();
      cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear> interpolator( *signedDistanceMap );
      
      // accumulate interpolated distances for this label
#pragma omp parallel for
      for ( int z = 0; z < referenceDims[2]; ++z )
	{
	cmtk::Vector3D v;
	cmtk::Types::DataItem dvalue;
	size_t i = z * referenceDims[0] * referenceDims[1];

	for ( int y = 0; y < referenceDims[1]; ++y )
	  {
	  for ( int x = 0; x < referenceDims[0]; ++x, ++i )
	    {
	    (*itX)->GetTransformedGrid( v, x, y, z );
	    if ( interpolator.GetDataAt( v, dvalue ) )
	      {
	      referenceInOutDistancePtr[i] += dvalue;
	      }
	    }
	  }
	}
      }

    // if this is not the first label, compare this label's sum distance map
    // (over all volumes) pixel by pixel and set this label where it is
    // closer than previous closest label
    if ( label )
      {
#ifdef CMTK_USE_GCD
      const cmtk::Threads::Stride stride( nPixelsReference );
      dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		      { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			  for ( int i = 0; i < static_cast<int>( nPixelsReference ); ++i )
#endif
	{
	if ( referenceInOutDistancePtr[i] < totalDistancePtr[i] )
	  {
	  totalDistancePtr[i] = referenceInOutDistancePtr[i];
	  resultPtr[i] = label;
	  }
	else
	  {
	  if ( !(referenceInOutDistancePtr[i] > totalDistancePtr[i]) )
	    {
	    resultPtr[i] = numLabels;
	    }
	  }
	}
#ifdef CMTK_USE_GCD
		      });
#endif

      }
    else
      {
      // for label 0, simply copy map.
// #pragma omp parallel for // not worth parallelizing a simple copy
      for ( size_t i = 0; i < nPixelsReference; ++i )
	{
	totalDistancePtr[i] = referenceInOutDistancePtr[i];
	}
      }

    if ( WriteDistanceMapNameMask )
      {
      cmtk::UniformVolume::SmartPtr distanceMap( (*volumes.begin())->CloneGrid() );
      distanceMap->SetData( totalDistance );
      
      char fname[PATH_MAX];
      if ( snprintf( fname, PATH_MAX, WriteDistanceMapNameMask, label ) > PATH_MAX )
	{
	cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
	}
      else
	{
	cmtk::VolumeIO::Write( *distanceMap, fname );
	}
      }

    if ( WriteLabelMapNameMask )
      {
      cmtk::UniformVolume::SmartPtr labelMap( (*volumes.begin())->CloneGrid() );
      labelMap->SetData( result );
      
      char fname[PATH_MAX];
      if ( snprintf( fname, PATH_MAX, WriteLabelMapNameMask, label ) > PATH_MAX )
	{
	cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
	}
      else
	{
	cmtk::VolumeIO::Write( *labelMap, fname );
	}
      }
    }

  return result;
}

cmtk::TypedArray::SmartPtr
AverageWindowed
( const cmtk::UniformVolume* referenceVolume,
  std::list<cmtk::UniformVolume::SmartPtr> volumes,
  std::list<cmtk::XformUniformVolume::SmartConstPtr> xforms,
  const short numLabels )
{
  const int distanceMapFlags = cmtk::UniformDistanceMap<float>::VALUE_WINDOW + cmtk::UniformDistanceMap<float>::SIGNED;

  const size_t nPixelsReference = referenceVolume->GetNumberOfPixels();
  const cmtk::DataGrid::IndexType& referenceDims = referenceVolume->GetDims();

  bool labelFlags[256];
  memset( labelFlags, 0, sizeof( labelFlags ) );

  for ( std::list<cmtk::UniformVolume::SmartPtr>::const_iterator itV = volumes.begin(); itV != volumes.end(); ++itV )
    {
    const cmtk::TypedArray* data = (*itV)->GetData();
#pragma omp parallel for
    for ( int i = 0; i < static_cast<int>( data->GetDataSize() ); ++i )
      {
      cmtk::Types::DataItem l;
      if ( data->Get( l, i ) )
	labelFlags[static_cast<short>( l )] = true;
      }
    }
  
  cmtk::TypedArray::SmartPtr result( cmtk::TypedArray::Create( cmtk::TYPE_FLOAT, nPixelsReference ) );
  result->BlockSet( 0 /*value*/, 0 /*idx*/, nPixelsReference /*len*/ );
  float* resultPtr = (float*) result->GetDataPtr();

  cmtk::TypedArray::SmartPtr resultDivider( cmtk::TypedArray::Create( cmtk::TYPE_SHORT, nPixelsReference ) );
  resultDivider->BlockSet( 0 /*value*/, 0 /*idx*/, nPixelsReference /*len*/ );
  short* resultDividerPtr = static_cast<short*>( resultDivider->GetDataPtr() );
  
  cmtk::TypedArray::SmartPtr countDistanceSamples( cmtk::TypedArray::Create( cmtk::TYPE_SHORT, nPixelsReference ) );
  short* countDistanceSamplesPtr = static_cast<short*>( countDistanceSamples->GetDataPtr() );
  
  cmtk::FloatArray::SmartPtr referenceInOutDistance( new cmtk::FloatArray( nPixelsReference ) );
  float* referenceInOutDistancePtr = referenceInOutDistance->GetDataPtrTemplate();

  cmtk::FloatArray::SmartPtr totalDistance( new cmtk::FloatArray( nPixelsReference ) );
  float* totalDistancePtr = totalDistance->GetDataPtrTemplate();

  totalDistance->BlockSet( 0 /*value*/, 0 /*idx*/, nPixelsReference /*len*/ );
  for ( int label = 0; label < numLabels; ++label )
    {
    /// skip labels that are not in any image.
    if ( ! labelFlags[label] ) continue;

    cmtk::DebugOutput( 1 ) << "Processing label #" << label << "\r";

    referenceInOutDistance->BlockSet( 0 /*value*/, 0 /*idx*/, nPixelsReference /*len*/ );
    countDistanceSamples->BlockSet( 0 /*value*/, 0 /*idx*/, nPixelsReference /*len*/ );

    std::list<cmtk::XformUniformVolume::SmartConstPtr>::const_iterator itX = xforms.begin();
    for ( std::list<cmtk::UniformVolume::SmartPtr>::const_iterator itV = volumes.begin(); itV != volumes.end(); ++itV, ++itX )
      {
      cmtk::UniformVolume::SmartPtr signedDistanceMap = cmtk::UniformDistanceMap<float>( *(*itV), distanceMapFlags, label, FeatureWindowRadius ).Get();
      cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear> interpolator( *signedDistanceMap );

      // accumulate interpolated distances for this label
#ifdef CMTK_USE_GCD
      const cmtk::Threads::Stride stride( referenceDims[2] );
      dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		      { for ( size_t z = stride.From( b ); z < stride.To( b ); ++z )
#else
#pragma omp parallel for
      for ( int z = 0; z < referenceDims[2]; ++z )
#endif
	{
	cmtk::Vector3D v;
	cmtk::Types::DataItem dvalue;

	size_t i = z * referenceDims[0] * referenceDims[1];
	for ( int y = 0; y < referenceDims[1]; ++y )
	  {
	  for ( int x = 0; x < referenceDims[0]; ++x, ++i )
	    {
	    (*itX)->GetTransformedGrid( v, x, y, z );
	    if ( interpolator.GetDataAt( v, dvalue ) )
	      {
	      referenceInOutDistancePtr[i] += dvalue;
	      ++countDistanceSamplesPtr[i];
	      }
	    }
	  }
	}
#ifdef CMTK_USE_GCD
		      });
#endif
      }
    
    // if this is not the first label, compare this label's sum distance map
    // (over all volumes) pixel by pixel and set this label where it is
    // closer than previous closest label
    if ( label )
      {
#ifdef CMTK_USE_GCD
      const cmtk::Threads::Stride stride( nPixelsReference );
      dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		      { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			  for ( int i = 0; i < static_cast<int>( nPixelsReference ); ++i )
#endif
	{
	if ( countDistanceSamplesPtr[i] )
	  {
	  referenceInOutDistancePtr[i] /= countDistanceSamplesPtr[i];
	  
	  if ( referenceInOutDistancePtr[i] < totalDistancePtr[i] )
	    {
	    totalDistancePtr[i] = referenceInOutDistancePtr[i];
	    resultPtr[i] = label;
	    resultDividerPtr[i] = 1;
	    }
	  else
	    {
	    if ( ! (referenceInOutDistancePtr[i] > totalDistancePtr[i] ) )
	      {
	      resultPtr[i] += label;
	      ++resultDividerPtr[i];
	      }
	    }
	  }
	}
#ifdef CMTK_USE_GCD
		      });
#endif
      }
    else
      {
      // for label 0, simply copy map.
#ifdef CMTK_USE_GCD
      const cmtk::Threads::Stride stride( nPixelsReference );
      dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		      { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			  for ( int i = 0; i < static_cast<int>( nPixelsReference ); ++i )
#endif
	{
	if ( countDistanceSamplesPtr[i] )
	  {
	  referenceInOutDistancePtr[i] /= countDistanceSamplesPtr[i];
	  totalDistancePtr[i] = referenceInOutDistancePtr[i];
	  resultDividerPtr[i] = 1;
	  }
	}
#ifdef CMTK_USE_GCD
		      });
#endif
      }
    }

  // compute average between min and max window result
#pragma omp parallel for
  for ( int i = 0; i < static_cast<int>( nPixelsReference ); ++i )
    {
    if ( resultDividerPtr[i] )
      {
      resultPtr[i] /= resultDividerPtr[i];
      }
    else
      {
      resultPtr[i] = numLabels;
      }
    resultDividerPtr[i] = 0;
    }

  result = cmtk::TypedArray::SmartPtr( result->Convert( cmtk::TYPE_SHORT ) );
  return result;
}

void
AddVolumeFile
( const char* fileName, std::list<cmtk::UniformVolume::SmartPtr>& volumeList )
{
  cmtk::DebugOutput( 1 ) << "Opening image " << fileName << ".\n";
  
  cmtk::UniformVolume::SmartPtr nextVolume( cmtk::VolumeIO::ReadOriented( fileName ) );
  
  if ( ! nextVolume || ! nextVolume->GetData() )
    {
    cmtk::StdErr << "WARNING: Could not read volume -- ignoring this.\n";
    }
  else
    {
    if ( PaddingFlag )
      {
      nextVolume->GetData()->SetPaddingValue( PaddingValue );
      }

    if ( ScaleIntensitiesToLabels )
      {
      nextVolume->GetData()->RescaleToRange( cmtk::Types::DataItemRange( 0, NumberOfLabels-1 ) );
      }
    
    if ( nextVolume->GetData()->GetType() != cmtk::TYPE_SHORT )
      {
      cmtk::StdErr << "WARNING: converting data to 'short'\n";
      
      nextVolume->SetData( cmtk::TypedArray::SmartPtr( nextVolume->GetData()->Convert( cmtk::TYPE_SHORT ) ) );
      }
    volumeList.push_back( nextVolume );
    }
}

void
AddVolumeStudyList
( const char* listName, std::list<cmtk::UniformVolume::SmartPtr>& volumeList,
  std::list<cmtk::XformUniformVolume::SmartConstPtr>& xformList,
  cmtk::UniformVolume::SmartPtr& referenceVolume )
{
  cmtk::DebugOutput( 1 ) << "Opening studylist " << listName << ".\n";

  cmtk::TypedStreamStudylist studylist;
  studylist.Read( listName );

  if ( ! referenceVolume )
    {
    const std::string actualPath = cmtk::StrReplace( studylist.GetReferenceStudyPath(), ReplaceMap );
    referenceVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadGridOriented( actualPath.c_str() ) );
    if ( ! referenceVolume )
      {
      cmtk::StdErr << "WARNING: could not read reference volume " << actualPath.c_str() << "\n";
      return;
      }
    }

  const std::string actualPath = cmtk::StrReplace( studylist.GetFloatingStudyPath(), ReplaceMap );
  cmtk::UniformVolume::SmartPtr floatingVolume( cmtk::VolumeIO::ReadOriented( actualPath.c_str() ) );

  if ( ! floatingVolume )
    {
    cmtk::StdErr << "WARNING: could not read floating volume " << actualPath.c_str() << " (skipping)\n";
    return;
    }

  if ( ! floatingVolume->GetData() )
    {
    cmtk::StdErr << "WARNING: could not read data from floating volume " << actualPath.c_str() << " (skipping)\n";
    return;
    }

  if ( PaddingFlag )
    {
    floatingVolume->GetData()->SetPaddingValue( PaddingValue );
    }

  if ( ScaleIntensitiesToLabels )
    {
    floatingVolume->GetData()->RescaleToRange( cmtk::Types::DataItemRange( 0, NumberOfLabels-1) );
    }

  if ( DownsampleVolumeStr )
    {
    floatingVolume = cmtk::UniformVolume::SmartPtr( floatingVolume->GetDownsampledAndAveraged( DownsampleVolume ) );
    cmtk::DebugOutput( 1 ) << "Downsampling atlas by factors " << DownsampleVolumeStr << "\n";
    }
  
  volumeList.push_back( floatingVolume );

  cmtk::SplineWarpXform::SmartConstPtr splinexform( cmtk::SplineWarpXform::SmartConstPtr::DynamicCastFrom( studylist.GetWarpXform() ) );
  if ( !splinexform )
    {
    cmtk::AffineXform::SmartConstPtr affinexform( studylist.GetAffineXform() );

    if ( !affinexform )
      {
      cmtk::StdErr << "WARNING: no transformation in studylist " << listName << "\n";
      return;
      }
    else
      {
      xformList.push_back( cmtk::XformUniformVolume::SmartConstPtr( new cmtk::AffineXformUniformVolume( *(referenceVolume), *(affinexform) ) ) );
      }
    }
  else
    {
    xformList.push_back( cmtk::XformUniformVolume::SmartConstPtr( new cmtk::SplineWarpXformUniformVolume( *(referenceVolume), splinexform ) ) );
    }
}

void
AddVolumeStudyList
( const char* listName, std::list<cmtk::UniformVolume::SmartPtr>& volumeList,
  cmtk::UniformVolume::SmartPtr& referenceVolume )
{
  cmtk::DebugOutput( 1 ) << "Opening studylist " << listName << ".\n";

  cmtk::TypedStreamStudylist studylist;
  studylist.Read( listName );

  if ( ! referenceVolume )
    {
    const std::string actualPath = cmtk::StrReplace( studylist.GetReferenceStudyPath(), ReplaceMap );
    referenceVolume = cmtk::UniformVolume::SmartPtr( cmtk::VolumeIO::ReadGridOriented( actualPath.c_str() ) );
    if ( ! referenceVolume )
      {
      cmtk::StdErr << "WARNING: could not read reference volume!";
      return;
      }
    }

  const std::string actualPath = cmtk::StrReplace( studylist.GetFloatingStudyPath(), ReplaceMap );
  cmtk::UniformVolume::SmartPtr floatingVolume( cmtk::VolumeIO::ReadOriented( actualPath.c_str() ) );

  if ( ! floatingVolume->GetData() )
    {
    cmtk::StdErr << "WARNING: could not read data from floating volume " << actualPath.c_str() << " (skipping)\n";
    return;
    }

  if ( PaddingFlag )
    {
    floatingVolume->GetData()->SetPaddingValue( PaddingValue );
    }

  if ( DownsampleVolumeStr )
    {
    floatingVolume = cmtk::UniformVolume::SmartPtr( floatingVolume->GetDownsampledAndAveraged( DownsampleVolume ) );
    cmtk::DebugOutput( 1 ) << "Downsampling atlas by factors " << DownsampleVolumeStr << "\n";
    }
  
  cmtk::ReformatVolume reformat;
  reformat.SetReferenceVolume( referenceVolume );
  reformat.SetFloatingVolume( floatingVolume );
  reformat.SetAffineXform( studylist.GetAffineXform() );
  reformat.SetWarpXform( studylist.GetWarpXform() );
  reformat.SetInterpolation( cmtk::Interpolators::PARTIALVOLUME );
  reformat.SetUsePaddingValue( true );
  reformat.SetPaddingValue( 255 );

  cmtk::UniformVolume::SmartPtr dataVolume( reformat.PlainReformat() );
  dataVolume->GetData()->ReplacePaddingData( 0 );

  if ( ScaleIntensitiesToLabels )
    {
    dataVolume->GetData()->RescaleToRange( cmtk::Types::DataItemRange( 0, NumberOfLabels-1 ) );
    }

  if ( dataVolume->GetData()->GetType() != cmtk::TYPE_SHORT )
    {
    cmtk::StdErr << "WARNING: converting data to 'short'\n";
    
    dataVolume->SetData( cmtk::TypedArray::SmartPtr( dataVolume->GetData()->Convert( cmtk::TYPE_SHORT ) ) );
    } 
  volumeList.push_back( dataVolume );
}

int
doMain ( const int argc, const char* argv[] ) 
{
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Shape-based averaging" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Average segmentations (label fields) using the Euclidean Distance Transform." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "average_edt [options] list0 [list1 ...]" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( "interpolate-image" ), &InterpolateDistance, false, "Interpolate target image, compute distance on reference grid" );
    cl.AddSwitch( Key( "interpolate-distance" ), &InterpolateDistance, true, "Compute distance maps on target image grid and interpolate distance" );

    cl.AddOption( Key( "feature-window-radius" ), &FeatureWindowRadius, "Radius of feature value window." );
    cl.AddSwitch( Key( "scale-to-labels" ), &ScaleIntensitiesToLabels, true, "Scale intensity images to given number of labels." );


    cl.AddOption( Key( "downsample-volume" ), &DownsampleVolumeStr, "Volume (atlas) downsampling factors 'dx,dy,dz'" );
    cl.AddOption( Key( 'n', "num-labels" ), &NumberOfLabels, "Number of labels. It is assumed that only values [0..num] occur in the images" );
    cl.AddOption( Key( 'p', "padding" ), &PaddingValue, "Padding value in input image", &PaddingFlag );

    cl.AddCallback( Key( "replace-from" ), AddReplaceFrom, "Replace from pattern" );
    cl.AddCallback( Key( "replace-to" ), AddReplaceTo, "Replace to pattern" );
    cl.AddOption( Key( 'o', "output" ), &OutputFileName, "File name for output segmentation file." );

    cl.AddOption( Key( "write-dmap" ), &WriteDistanceMapNameMask, "Write intermediate distance maps [file name mask]" );
    cl.AddOption( Key( "write-lmap" ), &WriteLabelMapNameMask, "Write intermediate label maps [file name mask]" );

    cl.Parse( argc, argv );

    const char* next = cl.GetNext();
    while ( next )
      {
      InputFileList.push_back( next );
      next = cl.GetNextOptional();
      }

    if ( DownsampleVolumeStr )
      {
      if ( 3 != sscanf( DownsampleVolumeStr, "%d,%d,%d", &DownsampleVolume[0], &DownsampleVolume[1], &DownsampleVolume[2] ) )
	{
	DownsampleVolume[0] = DownsampleVolume[1] = DownsampleVolume[2] = 1;
	DownsampleVolumeStr = NULL;
	}
      }
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }

  std::list<cmtk::UniformVolume::SmartPtr> volumeList;
  cmtk::UniformVolume::SmartPtr referenceVolume;
  cmtk::TypedArray::SmartPtr avgArray;

  if ( InterpolateDistance )
    {
    std::list<cmtk::XformUniformVolume::SmartConstPtr> xformList;
    for ( std::list<const char*>::const_iterator it = InputFileList.begin(); 
	  it != InputFileList.end(); ++it ) 
      {
      switch ( cmtk::FileFormat::Identify( *it ) )
	{
	case cmtk::FILEFORMAT_STUDYLIST:
	{
	AddVolumeStudyList( *it, volumeList, xformList, referenceVolume );
	break;
	}
	default:
	{
	cmtk::StdErr << "ERROR: all inputs must include transformation in distance map interpolation mode.\n";
	throw cmtk::ExitException( 1 );
	break;
	}
	}
      }

    const double timeBaseline = cmtk::Timers::GetTimeProcess();
    if ( FeatureWindowRadius > 0 )
      avgArray = cmtk::TypedArray::SmartPtr( AverageWindowed( referenceVolume, volumeList, xformList, NumberOfLabels ) );
    else
      avgArray = cmtk::TypedArray::SmartPtr( Average( referenceVolume, volumeList, xformList, NumberOfLabels ) );

    cmtk::DebugOutput( 1 ).GetStream().printf( "Time %f sec\n", cmtk::Timers::GetTimeProcess() - timeBaseline );
    }
  else
    {
    for ( std::list<const char*>::const_iterator it = InputFileList.begin(); it != InputFileList.end(); ++it ) 
      {
      switch ( cmtk::FileFormat::Identify( *it ) )
	{
	case cmtk::FILEFORMAT_STUDYLIST:
	{
	AddVolumeStudyList( *it, volumeList, referenceVolume );
	break;
	}
	default:
	{
	AddVolumeFile( *it, volumeList );
	break;
	}
	}
      }
    
    const double timeBaseline = cmtk::Timers::GetTimeProcess();
    if ( FeatureWindowRadius > 0 )
      avgArray = cmtk::TypedArray::SmartPtr( AverageWindowed( volumeList, NumberOfLabels ) );      
    else
      avgArray = cmtk::TypedArray::SmartPtr( Average( volumeList, NumberOfLabels ) );
    
    cmtk::DebugOutput( 1 ).GetStream().printf( "Time %f sec\n", cmtk::Timers::GetTimeProcess() - timeBaseline );
    }
    
  if ( OutputFileName )
    {
    cmtk::UniformVolume::SmartPtr volume = referenceVolume;
    if ( !volume )
      volume = *(volumeList.begin());

    volume->SetData( avgArray );
    cmtk::VolumeIO::Write( *volume, OutputFileName );
    }

  return 0;
}

#include "cmtkSafeMain"
