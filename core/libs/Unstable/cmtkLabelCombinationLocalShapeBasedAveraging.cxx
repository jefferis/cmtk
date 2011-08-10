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

#include "cmtkLabelCombinationLocalShapeBasedAveraging.h"

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
cmtk::LabelCombinationLocalShapeBasedAveraging::AddAtlas( const UniformVolume::SmartConstPtr image, const UniformVolume::SmartConstPtr atlas )
{
  Superclass::AddAtlas( image, atlas );

  this->m_AtlasDMaps.push_back( UniformDistanceMap<double>( *atlas, UniformDistanceMap<double>::SIGNED | UniformDistanceMap<double>::SQUARED ).Get() );
}

cmtk::TypedArray::SmartPtr 
cmtk::LabelCombinationLocalShapeBasedAveraging::GetResult() const
{
  const UniformVolume& targetImage = *(this->m_TargetImage);
  cmtk::TypedArray::SmartPtr result( TypedArray::Create( TYPE_SHORT, targetImage.GetNumberOfPixels() ) );
  
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
cmtk::LabelCombinationLocalShapeBasedAveraging::ComputeResultForRegion( const Self::TargetRegionType& region, TypedArray& result ) const
{
  const UniformVolume& targetImage = *(this->m_TargetImage);
  const size_t nAtlases = this->m_AtlasImages.size();
  std::vector<bool> valid( nAtlases );
  std::vector<short> labels( nAtlases );  
  std::vector<Types::DataItem> weights( nAtlases );  

  for ( RegionIndexIterator<TargetRegionType> it( region ); it != it.end(); ++it )
    {
    const size_t i = targetImage.GetOffsetFromIndex( it.Index() );

    for ( size_t n = 0; n < nAtlases; ++n )
      {
      Types::DataItem value;
      if ( (valid[n] = this->m_AtlasLabels[n]->GetData()->Get( value, i ) ) )
	labels[n] = static_cast<short>( value );
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
      // Compute weights for the atlases from local image patch similarity.
      const TargetRegionType patchRegion( Max( region.From(), it.Index() - this->m_PatchRadius ), Min( region.To(), it.Index() + this->m_PatchRadius ) );
      TypedArray::SmartConstPtr targetDataPatch( targetImage.GetRegionData( patchRegion ) );

      Types::DataItem minWeight = FLT_MAX;
      Types::DataItem maxWeight = FLT_MIN;
      for ( size_t n = 0; n < nAtlases; ++n )
	{
	if ( valid[n] )
	  {
	  TypedArray::SmartConstPtr atlasDataPatch( this->m_AtlasImages[n]->GetRegionData( patchRegion ) );
	  weights[n] = TypedArraySimilarity::GetCrossCorrelation( targetDataPatch, atlasDataPatch );
	  
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
	  totalDistance += (weights[n]-minWeight)/maxWeight * this->m_AtlasDMaps[n]->GetDataAt( i );
	  }
	}
	  
      result.Set( (totalDistance <= 0) ? 1 : 0, i );
      }
    }
}
