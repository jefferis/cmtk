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

#include "cmtkLabelCombinationLocalShapeBasedAveraging.h"

#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>
#include <System/cmtkDebugOutput.h>

#include <Base/cmtkRegionIndexIterator.h>
#include <Base/cmtkUniformDistanceMap.h>
#include <Base/cmtkMathFunctionWrappers.h>

#include <Registration/cmtkTypedArraySimilarity.h>

#ifdef _OPENMP
#  include <omp.h>
#endif

cmtk::TypedArray::SmartPtr 
cmtk::LabelCombinationLocalShapeBasedAveraging::GetResult() const
{
  const UniformVolume& targetImage = *(this->m_TargetImage);

  const size_t nPixels = targetImage.GetNumberOfPixels();
  cmtk::TypedArray::SmartPtr result( TypedArray::Create( TYPE_SHORT, nPixels ) );
  std::vector<float> resultDistance( nPixels, 1.0 );
  
  const TargetRegionType region = targetImage.CropRegion();
  
  // signed distance maps for the atlas label maps.
  const size_t nAtlases = this->m_AtlasImages.size();
  std::vector<UniformVolume::SmartConstPtr> atlasDMaps( nAtlases );

  const int maxLabelValue = (this->m_MaxLabelValue>0) ? this->m_MaxLabelValue : this->ComputeMaximumLabelValue();
  for ( int label = 0; label <= maxLabelValue; ++label )
    {
    if ( this->ComputeLabelNumberOfPixels( label ) > 0 ) // skip unused label values
      {
      DebugOutput( 2 ) << "Processing label " << label << "\n";

      DebugOutput( 5 ) << "  Creating distance maps\n";
      for ( size_t n = 0; n < nAtlases; ++n )
	{
	atlasDMaps[n] = ( UniformDistanceMap<float>( *(this->m_AtlasLabels[n]), DistanceMap::SIGNED | DistanceMap::SQUARED | DistanceMap::VALUE_EXACT, label ).Get() );
	}
      
      DebugOutput( 5 ) << "  Combining distance maps\n";
#ifdef _OPENMP
#pragma omp parallel for
      for ( int slice = region.From()[2]; slice < region.To()[2]; ++slice )
	{
	TargetRegionType threadRegion = region;
	threadRegion.From()[2] = slice;
	threadRegion.To()[2] = slice+1;
	
	this->ComputeResultForRegion( *result, resultDistance, label, threadRegion, atlasDMaps );
	}
#else // _OPENMP
      this->ComputeResultForRegion( *result, resultDistance, label, region, atlasDMaps );
#endif // _OPENMP
      }
    }
  
  return result;
}

void
cmtk::LabelCombinationLocalShapeBasedAveraging::ComputeResultForRegion
( TypedArray& result, std::vector<float>& resultDistance, const int label, const Self::TargetRegionType& region, std::vector<UniformVolume::SmartConstPtr> dmaps ) const
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
      if ( (valid[n] = dmaps[n]->GetData()->Get( value, i ) ) )
	labels[n] = static_cast<short>( (value <= 0) ? label : -1 );
      }
    
    // detect local outliers in the distance maps, ie., grossly misregistered atlases
    if ( this->m_DetectLocalOutliers )
      {
      // create vector of distance values
      size_t nn = 0;
      for ( size_t n = 0; n < nAtlases; ++n )
	{
	if ( valid[n] )
	  distances[nn++] = dmaps[n]->GetDataAt( i );
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
	  const double d = dmaps[n]->GetDataAt( i );
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
	totalDistance += (weights[n]-minWeight)/maxWeight * dmaps[n]->GetDataAt( i + bestPatchOffset[n] );
	}
      }  
    
    if ( totalDistance < resultDistance[i] )
      {
      result.Set( label, i );
      resultDistance[i] = static_cast<float>( totalDistance );
      }
    }
}
