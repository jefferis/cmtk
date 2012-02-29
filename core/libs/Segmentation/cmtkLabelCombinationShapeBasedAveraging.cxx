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

#include "cmtkLabelCombinationShapeBasedAveraging.h"

#include <System/cmtkDebugOutput.h>
#include <System/cmtkThreads.h>

#include <Base/cmtkDistanceMap.h>
#include <Base/cmtkUniformDistanceMap.h>
#include <Base/cmtkTypedArray.h>
#include <Base/cmtkTemplateArray.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

LabelCombinationShapeBasedAveraging::LabelCombinationShapeBasedAveraging( const std::vector<UniformVolume::SmartConstPtr>& labelImages, const Self::LabelIndexType numberOfLabels )
  : m_NumberOfLabels( numberOfLabels ),
    m_LabelImages( labelImages )
{
  if ( ! this->m_NumberOfLabels )
    {
    this->m_NumberOfLabels = 1;
    for ( size_t k = 0; k < this->m_LabelImages.size(); ++k )
      {
      const Types::DataItemRange range = this->m_LabelImages[k]->GetData()->GetRange();
      this->m_NumberOfLabels = std::max( this->m_NumberOfLabels, static_cast<Self::LabelIndexType>( 1 + range.m_UpperBound ) );
      }
    
    DebugOutput( 9 ) << "Determined number of labels to be " << this->m_NumberOfLabels << "\n";
    }

  this->m_NumberOfPixels = this->m_LabelImages[0]->GetNumberOfPixels();
  
  this->m_LabelFlags.resize( this->m_NumberOfLabels, false );
  for ( size_t k = 0; k < this->m_LabelImages.size(); ++k )
    {
    const cmtk::TypedArray& data = *(this->m_LabelImages[k]->GetData());
    
    cmtk::Types::DataItem l;
    for ( size_t i = 0; i < this->m_NumberOfPixels; ++i )
      {
      if ( data.Get( l, i ) )
	this->m_LabelFlags[static_cast<unsigned short>( l )] = true;
      }
    }
  
}

TypedArray::SmartPtr
LabelCombinationShapeBasedAveraging::GetResult( const bool detectOutliers ) const
{
  cmtk::TypedArray::SmartPtr result( cmtk::TypedArray::Create( cmtk::TYPE_USHORT, this->m_NumberOfPixels ) );
  result->BlockSet( 0 /*value*/, 0 /*idx*/, this->m_NumberOfPixels /*len*/ );
  
  std::vector<Self::DistanceMapRealType> totalDistance( this->m_NumberOfPixels, 0.0 );
  std::vector<Self::DistanceMapRealType> labelDistanceMap( this->m_NumberOfPixels );

  for ( int label = 0; label < this->m_NumberOfLabels; ++label )
    {
    /// skip labels that are not in any image.
    if ( ! this->m_LabelFlags[label] ) continue;

    cmtk::DebugOutput( 1 ) << "Processing label #" << label << "\r";

    std::fill( labelDistanceMap.begin(), labelDistanceMap.end(), 0.0 );

    if ( detectOutliers )
      {
      // if this is the first label, write directly to totalDistance, otherwise labelDistanceMap
      this->ProcessLabelExcludeOutliers( label, (label == 0) ? totalDistance : labelDistanceMap );
      }
    else
      {
      // if this is the first label, write directly to totalDistance, otherwise labelDistanceMap
      this->ProcessLabelIncludeOutliers( label, (label == 0) ? totalDistance : labelDistanceMap );
      }
    
    // compare sum over all inputs of this label's distance maps pixel by pixel with current total distance map.
    // Set result map to this label where it is closer than previous closest label.
    if ( label )
      {
#pragma omp parallel for
      for ( int i = 0; i < static_cast<int>( this->m_NumberOfPixels ); ++i )
	{
	if ( labelDistanceMap[i] < totalDistance[i] )
	  {
	  totalDistance[i] = labelDistanceMap[i];
	  result->Set( label, i );
	  }
	else
	  {
	  if ( !(labelDistanceMap[i] > totalDistance[i]) )
	    {
	    result->Set( this->m_NumberOfLabels, i );
	    }	  
	  }
	}
      }
    }
  
  return result;
}

void
LabelCombinationShapeBasedAveraging::ProcessLabelExcludeOutliers
( const Self::LabelIndexType label, std::vector<Self::DistanceMapRealType>& labelDistanceMap ) const
{
  const size_t nLabelMaps = this->m_LabelImages.size();
  std::vector<cmtk::UniformVolume::SmartConstPtr> signedDistanceMaps( nLabelMaps );
  for ( size_t k = 0; k < nLabelMaps; ++k )
    {
    signedDistanceMaps[k] = cmtk::UniformDistanceMap<Self::DistanceMapRealType>( *(this->m_LabelImages[k]), DistanceMap::VALUE_EXACT + DistanceMap::SIGNED, label ).Get();
    }

  std::vector<Self::DistanceMapRealType> distances( nLabelMaps );
  
  for ( int i = 0; i < static_cast<int>( this->m_NumberOfPixels ); ++i )
    {
    for ( size_t k = 0; k < nLabelMaps; ++k )
      {
      distances[k] = signedDistanceMaps[k]->GetDataAt( i );
      }
    
    // sort distance
    std::sort( distances.begin(), distances.end() );
    
    // determine 1st and 3rd quartile values
    const double Q1 = distances[static_cast<size_t>( 0.25 * nLabelMaps )];
    const double Q3 = distances[static_cast<size_t>( 0.75 * nLabelMaps )];
    
    // compute thresholds from quartiles and inter-quartile range
    const double lThresh = Q1 - 1.5 * (Q3-Q1);
    const double uThresh = Q3 + 1.5 * (Q3-Q1);
    
    for ( size_t k = 0; k < nLabelMaps; ++k )
      {
      if ( (distances[k] >= lThresh) && (distances[k] <= uThresh) )
	labelDistanceMap[i] += distances[k];
      }
    }
}

void
LabelCombinationShapeBasedAveraging::ProcessLabelIncludeOutliers
( const Self::LabelIndexType label, std::vector<Self::DistanceMapRealType>& labelDistanceMap ) const
{
  for ( size_t k = 0; k < this->m_LabelImages.size(); ++k )
    {
    cmtk::UniformVolume::SmartPtr signedDistanceMap = cmtk::UniformDistanceMap<Self::DistanceMapRealType>( *(this->m_LabelImages[k]), DistanceMap::VALUE_EXACT + DistanceMap::SIGNED, label ).Get();
    const Self::DistanceMapRealType* signedDistancePtr = static_cast<const Self::DistanceMapRealType*>( signedDistanceMap->GetData()->GetDataPtr() );
    
#pragma omp parallel for
    for ( int i = 0; i < static_cast<int>( this->m_NumberOfPixels ); ++i )
      {
      labelDistanceMap[i] += signedDistancePtr[i];
      }
    }
}

} // namespace cmtk
