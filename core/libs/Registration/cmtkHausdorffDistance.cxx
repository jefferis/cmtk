/*
//
//  Copyright 2011 SRI International
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

#include "cmtkHausdorffDistance.h"

#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkUniformDistanceMap.h>

#include <algorithm>

cmtk::HausdorffDistance::HausdorffDistance( UniformVolume::SmartConstPtr& image0, UniformVolume::SmartConstPtr& image1 )
  : m_Image0( image0 ),
    m_Image1( image1 )
{
  if ( !this->m_Image0->GridMatches( *(this->m_Image1) ) )
    {
    cmtk::StdErr << "ERROR: the two image grids don't match.\n";
    throw cmtk::ExitException( 1 );
    }
}

cmtk::Types::Coordinate 
cmtk::HausdorffDistance::GetBinary() const 
{
  typedef UniformDistanceMap<Types::Coordinate> DistanceMapType;

  UniformVolume::SmartConstPtr distance0 = DistanceMapType( *(this->m_Image0), DistanceMapType::DEFAULT ).Get(); 
  UniformVolume::SmartConstPtr distance1 = DistanceMapType( *(this->m_Image1), DistanceMapType::DEFAULT ).Get(); 

  return std::max( Self::HalfDistanceBinary( *(this->m_Image0), *distance1 ), Self::HalfDistanceBinary( *(this->m_Image1), *distance0 ) );
}

cmtk::Types::Coordinate 
cmtk::HausdorffDistance::HalfDistanceBinary( const UniformVolume& image, const UniformVolume& dmap )
{
  Types::Coordinate maxDistance = 0;

  const size_t nPixels = image.GetNumberOfPixels();
  for ( size_t n = 0; n < nPixels; ++n )
    {
    if ( image.GetDataAt( n ) )
      {
      maxDistance = std::max( maxDistance, dmap.GetDataAt( n ) );
      }
    }

  return maxDistance;
}
