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

#include <IO/cmtkVolumeIO.h>

#include <Registration/cmtkReformatVolume.h>

#include <math.h>
#include <vector>

#ifdef CMTK_USE_GCD
#  include <dispatch/dispatch.h>
#endif

const char *OutputFileName = "sba.nii";

std::vector<std::string> InputFileVector;

unsigned short NumberOfLabels = 0;

bool PaddingFlag = false;
float PaddingValue = 0;

cmtk::TypedArray::SmartPtr
Average
( std::vector<cmtk::UniformVolume::SmartPtr> volumes, const unsigned short numLabels )
{
  const int distanceMapFlags = cmtk::UniformDistanceMap<float>::VALUE_EXACT + cmtk::UniformDistanceMap<float>::SIGNED;

  const size_t numPixels = (*volumes.begin())->GetNumberOfPixels();

  std::vector<bool> labelFlags( numLabels, false );  
  for ( std::vector<cmtk::UniformVolume::SmartPtr>::const_iterator it = volumes.begin(); it != volumes.end(); ++it )
    {
    const cmtk::TypedArray& data = *((*it)->GetData());

    cmtk::Types::DataItem l;
    for ( size_t i = 0; i < numPixels; ++i )
      {
      if ( data.Get( l, i ) )
	labelFlags[static_cast<unsigned short>( l )] = true;
      }
    }

  cmtk::TypedArray::SmartPtr result( cmtk::TypedArray::Create( cmtk::TYPE_USHORT, numPixels ) );
  result->BlockSet( 0 /*value*/, 0 /*idx*/, numPixels /*len*/ );
  unsigned short* resultPtr = static_cast<unsigned short*>( result->GetDataPtr() );
  
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

    for ( std::vector<cmtk::UniformVolume::SmartPtr>::const_iterator it = volumes.begin(); it != volumes.end(); ++it )
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
    }
  
  return result;
}

void
AddVolumeFile
( const char* fileName, std::vector<cmtk::UniformVolume::SmartPtr>& volumeVector )
{
  cmtk::DebugOutput( 1 ) << "Opening image " << fileName << ".\n";
  
  cmtk::UniformVolume::SmartPtr nextVolume( cmtk::VolumeIO::ReadOriented( fileName ) );
  
  if ( PaddingFlag )
    {
    nextVolume->GetData()->SetPaddingValue( PaddingValue );
    }
  
  if ( nextVolume->GetData()->GetType() != cmtk::TYPE_USHORT )
    {
    cmtk::StdErr << "WARNING: converting data to 'unsigned short'\n";
    
    nextVolume->SetData( cmtk::TypedArray::SmartPtr( nextVolume->GetData()->Convert( cmtk::TYPE_USHORT ) ) );
    }
  volumeVector.push_back( nextVolume );
}

int
doMain ( const int argc, const char* argv[] ) 
{
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Shape-based Averaging of label images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Average segmentations (label fields) using the Euclidean Distance Transform. All input images must be in the same space. EDT is computed in this space also. "
		       "See http://dx.doi.org/10.1109/TIP.2006.884936 for details of the underlying algorithm." );

    cl.AddParameterVector( &InputFileVector, "InputImageVector", "Input image file names." );

    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Input", "Input Options" );
    cl.AddOption( Key( 'n', "num-labels" ), &NumberOfLabels, "Number of labels. It is assumed that only values [0..num] occur in the images" );
    cl.AddOption( Key( 'p', "padding" ), &PaddingValue, "Padding value in input image", &PaddingFlag );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( 'o', "output" ), &OutputFileName, "File name for output segmentation file." );
    cl.EndGroup();

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }

  std::vector<cmtk::UniformVolume::SmartPtr> volumeVector;
  for ( std::vector<std::string>::const_iterator it = InputFileVector.begin(); it != InputFileVector.end(); ++it ) 
    {
    AddVolumeFile( it->c_str(), volumeVector );
    }

  const double timeBaseline = cmtk::Timers::GetTimeProcess();
  cmtk::TypedArray::SmartPtr avgArray = cmtk::TypedArray::SmartPtr( Average( volumeVector, NumberOfLabels ) );
  cmtk::DebugOutput( 1 ).GetStream().printf( "Time %f sec\n", cmtk::Timers::GetTimeProcess() - timeBaseline );
    
  cmtk::UniformVolume& volume = *(volumeVector[0]);
  volume.SetData( avgArray );
  cmtk::VolumeIO::Write( volume, OutputFileName );

  return 0;
}

#include "cmtkSafeMain"
