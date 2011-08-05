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

#include "cmtkLabelCombinationLocalVoting.h"

#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkRegionIndexIterator.h>

#include <Registration/cmtkTypedArraySimilarity.h>

void
cmtk::LabelCombinationLocalVoting::AddAtlas
( const UniformVolume::SmartConstPtr image, const UniformVolume::SmartConstPtr atlas )
{
  if ( !this->m_TargetImage->GridMatches( *image ) )
    {
    StdErr << "Atlas intensity image grid does not match target image.\n";
    throw ExitException( 1 );
    }

  this->m_AtlasImages.push_back( image );

  if ( !this->m_TargetImage->GridMatches( *atlas ) )
    {
    StdErr << "Atlas label image grid does not match target image.\n";
    throw ExitException( 1 );
    }

  this->m_AtlasLabels.push_back( atlas );
}

cmtk::TypedArray::SmartPtr 
cmtk::LabelCombinationLocalVoting::GetResult() const
{
  const UniformVolume& targetImage = *(this->m_TargetImage);
  cmtk::TypedArray::SmartPtr result( TypedArray::Create( TYPE_SHORT, targetImage.GetNumberOfPixels() ) );
  
  const size_t nAtlases = this->m_AtlasImages.size();

  std::vector<bool> valid( nAtlases );
  std::vector<short> labels( nAtlases );  
  std::vector<Types::DataItem> weights( nAtlases );

  const UniformVolume::RegionType region = targetImage.CropRegion();
  for ( RegionIndexIterator<UniformVolume::RegionType> it( region ); it != region.end(); ++it )
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
      result->SetPaddingAt( i );
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
      result->Set( labels[firstValid], i );
      }
    else
      {
      // Compute weights for the atlases from local image patch similarity.
      
      }
    }

  return result;
}