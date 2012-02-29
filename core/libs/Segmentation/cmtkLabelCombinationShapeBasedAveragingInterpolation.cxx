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

#include "cmtkLabelCombinationShapeBasedAveragingInterpolation.h"

#include <System/cmtkDebugOutput.h>
#include <System/cmtkThreads.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkDistanceMap.h>
#include <Base/cmtkUniformDistanceMap.h>
#include <Base/cmtkTypedArray.h>
#include <Base/cmtkTemplateArray.h>
#include <Base/cmtkLinearInterpolator.h>
#include <Base/cmtkUniformVolumeInterpolator.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

LabelCombinationShapeBasedAveragingInterpolation::LabelCombinationShapeBasedAveragingInterpolation
( const std::vector<UniformVolume::SmartConstPtr>& labelImages, const std::vector<cmtk::XformUniformVolume::SmartConstPtr> xformsToLabelImages, const UniformVolume::SmartConstPtr& targetGrid, 
  const Self::LabelIndexType numberOfLabels )
  : LabelCombinationShapeBasedAveraging( labelImages, numberOfLabels ),
    m_TargetGrid( targetGrid ),
    m_Transformations( xformsToLabelImages )
{
  if ( this->m_LabelImages.size() != this->m_Transformations.size() )
    {
    StdErr << "ERROR: number of transformations does not match number of input images\n";
    throw ExitException( 1 );
    }

  this->m_NumberOfPixels = this->m_TargetGrid->GetNumberOfPixels();
}

TypedArray::SmartPtr
LabelCombinationShapeBasedAveragingInterpolation::GetResult( const bool detectOutliers ) const
{
  const cmtk::DataGrid::IndexType& targetDims = this->m_TargetGrid->GetDims();

  cmtk::TypedArray::SmartPtr result( cmtk::TypedArray::Create( cmtk::TYPE_USHORT,this->m_NumberOfPixels ) );
  result->BlockSet( 0 /*value*/, 0 /*idx*/,this->m_NumberOfPixels /*len*/ );
  unsigned short* resultPtr = static_cast<unsigned short*>( result->GetDataPtr() );
  
  cmtk::FloatArray::SmartPtr referenceInOutDistance( new cmtk::FloatArray(this->m_NumberOfPixels ) );
  float* referenceInOutDistancePtr = referenceInOutDistance->GetDataPtrTemplate();

  cmtk::FloatArray::SmartPtr totalDistance ( new cmtk::FloatArray(this->m_NumberOfPixels ) );
  float* totalDistancePtr = totalDistance->GetDataPtrTemplate();

  totalDistance->BlockSet( 0 /*value*/, 0 /*idx*/,this->m_NumberOfPixels /*len*/ );
  for ( int label = 0; label < this->m_NumberOfLabels; ++label )
    {
    /// skip labels that are not in any image.
    if ( ! this->m_LabelFlags[label] ) continue;

    cmtk::DebugOutput( 1 ) << "Processing label #" << label << "\r";

    referenceInOutDistance->BlockSet( 0 /*value*/, 0 /*idx*/,this->m_NumberOfPixels /*len*/ );

    for ( size_t k = 0; k < this->m_LabelImages.size(); ++k )
      {
      cmtk::UniformVolume::SmartPtr signedDistanceMap = cmtk::UniformDistanceMap<float>( *(this->m_LabelImages[k]), cmtk::DistanceMap::VALUE_EXACT + cmtk::DistanceMap::SIGNED, label ).Get();
      cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear> interpolator( *signedDistanceMap );
      
      // accumulate interpolated distances for this label
#pragma omp parallel for
      for ( int z = 0; z < targetDims[2]; ++z )
	{
	cmtk::Vector3D v;
	cmtk::Types::DataItem dvalue;
	size_t i = z * targetDims[0] * targetDims[1];

	for ( int y = 0; y < targetDims[1]; ++y )
	  {
	  for ( int x = 0; x < targetDims[0]; ++x, ++i )
	    {
	    this->m_Transformations[k]->GetTransformedGrid( v, x, y, z );
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
      const cmtk::Threads::Stride stride(this->m_NumberOfPixels );
      dispatch_apply( stride.NBlocks(), dispatch_get_global_queue(0, 0), ^(size_t b)
		      { for ( size_t i = stride.From( b ); i < stride.To( b ); ++i )
#else
#pragma omp parallel for
			  for ( int i = 0; i < static_cast<int>(this->m_NumberOfPixels ); ++i )
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
	    resultPtr[i] = this->m_NumberOfLabels;
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
      for ( size_t i = 0; i <this->m_NumberOfPixels; ++i )
	{
	totalDistancePtr[i] = referenceInOutDistancePtr[i];
	}
      }
    }
  
  return result;
}

void
LabelCombinationShapeBasedAveragingInterpolation::ProcessLabelExcludeOutliers
( const Self::LabelIndexType label, std::vector<Self::DistanceMapRealType>& labelDistanceMap ) const
{
}

void
LabelCombinationShapeBasedAveragingInterpolation::ProcessLabelIncludeOutliers
( const Self::LabelIndexType label, std::vector<Self::DistanceMapRealType>& labelDistanceMap ) const
{
}

} // namespace cmtk
