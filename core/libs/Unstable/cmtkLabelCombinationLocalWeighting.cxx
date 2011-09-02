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

#include "cmtkLabelCombinationLocalWeighting.h"

#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkRegionIndexIterator.h>

#include <Registration/cmtkTypedArraySimilarity.h>

#include <vector>
#include <algorithm>

#ifdef _OPENMP
#  include <omp.h>
#endif

void
cmtk::LabelCombinationLocalWeighting::AddAtlasImage
( const UniformVolume::SmartConstPtr image )
{
  if ( !this->m_TargetImage->GridMatches( *image ) )
    {
    StdErr << "Atlas intensity image grid does not match target image.\n";
    throw ExitException( 1 );
    }

  this->m_AtlasImages.push_back( image );
}

void
cmtk::LabelCombinationLocalWeighting::ExcludeGlobalOutliers()
{
  std::vector<Types::DataItem> ncc( this->m_AtlasImages.size() );

  for ( size_t i = 0; i < this->m_AtlasImages.size(); ++i )
    {
    ncc[i] = TypedArraySimilarity::GetCrossCorrelation( this->m_TargetImage->GetData(), this->m_AtlasImages[i]->GetData() );
    }

  std::vector<Types::DataItem> nccSort = ncc;
  std::sort( nccSort.begin(), nccSort.end() );

  // determine 1st and 3rd quartile values
  const Types::DataItem Q1 = nccSort[static_cast<size_t>( 0.25 * nccSort.size() )];
  const Types::DataItem Q3 = nccSort[static_cast<size_t>( 0.75 * nccSort.size() )];
  
  // compute threshold from 1st quartile and inter-quartile range
  const Types::DataItem lThresh = Q1 - 1.5 * (Q3-Q1);

  size_t iAtlas = 0;
  for ( size_t i = 0; i < this->m_AtlasImages.size(); ++i )
    {
    if ( ncc[i] < lThresh )
      {
      this->DeleteAtlas( iAtlas );      
      }
    else
      {
      ++iAtlas;
      }
    }
}
