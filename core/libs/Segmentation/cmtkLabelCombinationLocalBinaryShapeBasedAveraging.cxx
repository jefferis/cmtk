/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012, 2014 SRI International
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

#include "cmtkLabelCombinationLocalBinaryShapeBasedAveraging.h"

#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkRegionIndexIterator.h>
#include <Base/cmtkUniformDistanceMap.h>
#include <Base/cmtkMathFunctionWrappers.h>

#include <Registration/cmtkTypedArraySimilarity.h>

#ifdef _OPENMP
#  include <omp.h>
#endif

void
cmtk::LabelCombinationLocalBinaryShapeBasedAveraging::AddAtlas( const UniformVolume::SmartConstPtr image, const UniformVolume::SmartConstPtr atlas )
{
  Superclass::AddAtlasImage( image );

  this->m_AtlasDMaps.push_back( UniformDistanceMap<float>( *atlas, UniformDistanceMap<float>::SIGNED | UniformDistanceMap<float>::SQUARED ).Get() );
}

cmtk::TypedArray::SmartPtr 
cmtk::LabelCombinationLocalBinaryShapeBasedAveraging::GetResult() const
{
  const UniformVolume& targetImage = *(this->m_TargetImage);
  cmtk::TypedArray::SmartPtr result( TypedArray::Create( TYPE_SHORT, targetImage.GetNumberOfPixels() ) );
  result->SetDataClass( DATACLASS_LABEL );
  
  const TargetRegionType region = targetImage.CropRegion();

#ifdef _OPENMP
#pragma omp parallel for
  for ( int slice = region.From()[2]; slice < region.To()[2]; ++slice )
    {
    TargetRegionType threadRegion = region;
    threadRegion.From()[2] = slice;
    threadRegion.To()[2] = slice+1;
    
    this->ComputeResultForRegion( threadRegion, *result );
    }
#else // _OPENMP
  this->ComputeResultForRegion( region, *result );
#endif // _OPENMP
  
  return result;
}

void
cmtk::LabelCombinationLocalBinaryShapeBasedAveraging::ComputeResultForRegion( const Self::TargetRegionType& region, TypedArray& result ) const
{
  const UniformVolume& targetImage = *(this->m_TargetImage);
  const Self::TargetRegionType wholeImageRegion = targetImage.CropRegion();
  const size_t nAtlases = this->m_AtlasImages.size();
  std::vector<bool> valid( nAtlases );
  std::vector<short> labels( nAtlases );  
  std::vector<Types::DataItem> weights( nAtlases );  
  std::vector<size_t> bestPatchOffset( nAtlases );
  std::vector<double> distances( nAtlases );

  for ( RegionIndexIterator<TargetRegionType> it( region ); it != it.end(); ++it )
    {
    const size_t i = targetImage.GetOffsetFromIndex( it.Index() );

    for ( size_t n = 0; n < nAtlases; ++n )
      {
      Types::DataItem value;
      if ( (valid[n] = this->m_AtlasDMaps[n]->GetData()->Get( value, i ) ) )
	labels[n] = static_cast<short>( (value <= 0) ? 1 : 0 );
      }

    // detect local outliers in the distance maps, ie., grossly misregistered atlases
    if ( this->m_DetectLocalOutliers )
      {
      // create vector of distance values
      size_t nn = 0;
      for ( size_t n = 0; n < nAtlases; ++n )
	{
	if ( valid[n] )
	  distances[nn++] = this->m_AtlasDMaps[n]->GetDataAt( i );	
	}
      
      // sort distance
      std::sort( distances.begin(), distances.begin()+nn );

      // determine 1st and 3rd quartile values
      const double Q1 = distances[static_cast<size_t>( 0.25 * nn )];
      const double Q3 = distances[static_cast<size_t>( 0.75 * nn )];
      
      // compute thresholds from quartiles and inter-quartile range
      const double lThresh = Q1 - 1.5 * (Q3-Q1);
      const double uThresh = Q3 + 1.5 * (Q3-Q1);
      
      // mark as invalid those atlases with values outside the "inlier" range
      for ( size_t n = 0; n < nAtlases; ++n )
	{
	if ( valid[n] )
	  {
	  const double d = this->m_AtlasDMaps[n]->GetDataAt( i );
	  if ( (d < lThresh) || (d > uThresh) )
	    valid[n] = false;
	  }
	}      
      }
    
    // find first non-padding atlas label
    size_t firstValid = 0;
    while ( (firstValid < nAtlases) && !valid[firstValid] )
      ++firstValid;
    
    // if all input atlases are undefined (padding) for this pixel, set output to padding and skip to next pixel.
    if ( firstValid == nAtlases )
      {
      result.SetPaddingAt( i );
      continue;
      }

    // check if all (valid) input atlas labels are the same
    bool allTheSame = true;
    for ( size_t n = 1; n < nAtlases; ++n )
      {
      if ( valid[n] )
	{
	if ( labels[n] != labels[firstValid] )
	  {
	  allTheSame = false;
	  break;
	  }
	}
      }
    
    // no need for weighted combination if all labels are the same.
    if ( allTheSame )
      {
      result.Set( labels[firstValid] ? 1 : 0, i );
      }
    else
      {
      std::fill( weights.begin(), weights.end(), -1 );
      std::fill( bestPatchOffset.begin(), bestPatchOffset.end(), 0 );

      const TargetRegionType patchSearchRegion( Max( (-1)*wholeImageRegion.From(), this->m_SearchRegion.From() ), Min( wholeImageRegion.To() - it.Index(), this->m_SearchRegion.To() ) );
      for ( RegionIndexIterator<TargetRegionType> searchIt( patchSearchRegion ); searchIt != searchIt.end(); ++searchIt )
	{
	const TargetRegionType patchRegion( Max( wholeImageRegion.From(), it.Index() + searchIt.Index() - this->m_PatchRadius ), Min( wholeImageRegion.To(), it.Index() + searchIt.Index() + this->m_PatchRadiusPlusOne ) );
	TypedArray::SmartConstPtr targetDataPatch( targetImage.GetRegionData( patchRegion ) );
	
	for ( size_t n = 0; n < nAtlases; ++n )
	  {
	  if ( valid[n] )
	    {
	    TypedArray::SmartConstPtr atlasDataPatch( this->m_AtlasImages[n]->GetRegionData( patchRegion ) );
	    const Types::DataItem w = TypedArraySimilarity::GetCrossCorrelation( targetDataPatch, atlasDataPatch );
	    if ( w > weights[n] )
	      {
	      weights[n] = w;
	      bestPatchOffset[n] = targetImage.GetOffsetFromIndex( searchIt.Index() );
	      }
	    }
	  }
	}
      
      // Compute weights for the atlases from local image patch similarity.
      Types::DataItem minWeight = FLT_MAX;
      Types::DataItem maxWeight = FLT_MIN;
      
      for ( size_t n = 0; n < nAtlases; ++n )
	{
	if ( valid[n] )
	  {
	  minWeight = std::min( minWeight, weights[n] );
	  maxWeight = std::max( maxWeight, weights[n] );
	  }
	}
	
      maxWeight -= minWeight; // turn "max" into "range"
      
      double totalDistance = 0;
      for ( size_t n = 0; n < nAtlases; ++n )
	{
	if ( valid[n] )
	  {
	  totalDistance += (weights[n]-minWeight)/maxWeight * this->m_AtlasDMaps[n]->GetDataAt( i + bestPatchOffset[n] );
	  }
	}  
      
      result.Set( (totalDistance <= 0) ? 1 : 0, i );
      }
    }
}
