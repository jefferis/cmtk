/*
//
//  Copyright 2012 SRI International
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

#include "cmtkDetectPhantomMagphanEMR051.h"

#include <Base/cmtkMagphanEMR051.h>
#include <Base/cmtkUniformVolumePainter.h>

cmtk::DetectPhantomMagphanEMR051::DetectPhantomMagphanEMR051( UniformVolume::SmartConstPtr& phantomImage ) 
  : m_PhantomImage( phantomImage )
#ifdef CMTK_USE_FFTW
  , m_SphereDetector( *phantomImage )
#endif
{
}

std::vector<cmtk::UniformVolume::SpaceVectorType> 
cmtk::DetectPhantomMagphanEMR051::GetLandmarks()
{
  std::vector<cmtk::UniformVolume::SpaceVectorType> landmarks( MagphanEMR051::NumberOfSpheres );

  UniformVolume::SmartPtr excludeMask( this->m_PhantomImage->CloneGrid() );
  excludeMask->CreateDataArray( TYPE_BYTE, true /*setToZero*/ );

  // Find 60mm SNR sphere
  this->FindSpheres( landmarks.begin(), 1 /*n=1*/, MagphanEMR051::SphereTable[0].m_Diameter / 2, excludeMask );

  // Find the four 30mm spheres
  this->FindSpheres( landmarks.begin()+1, 4 /*n=4*/, MagphanEMR051::SphereTable[1].m_Diameter / 2, excludeMask );

  // Find the two 15mm spheres
  this->FindSpheres( landmarks.begin()+5, 2 /*n=2*/,  MagphanEMR051::SphereTable[5].m_Diameter / 2, excludeMask );

  return landmarks;
}

void
cmtk::DetectPhantomMagphanEMR051::FindSpheres( std::vector<cmtk::UniformVolume::SpaceVectorType>::iterator dest, const int nSpheres, const Types::Coordinate radius, UniformVolume::SmartPtr& excludeMask )
{
  UniformVolumePainter maskPainter( excludeMask, UniformVolumePainter::COORDINATES_ABSOLUTE );

  TypedArray::SmartPtr filterResponse( this->m_SphereDetector.GetFilteredImageData( radius ) );
  for ( int n = 0; n < nSpheres; ++n )
    {
    size_t maxIndex = 0;
    Types::DataItem maxValue = filterResponse->ValueAt( 0 );

    for ( size_t px = 0; px < filterResponse->GetDataSize(); ++px )
      {
      if ( excludeMask->GetDataAt( px ) == 0 )
	{
	const Types::DataItem value = filterResponse->ValueAt( px );
	if ( value > maxValue )
	  {
	  maxValue = value;
	  maxIndex = px;
	  }
	}
      }
    
    cmtk::UniformVolume::SpaceVectorType location = this->m_PhantomImage->GetGridLocation( maxIndex );
    
    // update exclusion mask
    maskPainter.DrawSphere( location, radius, 1 );
    
    // put landmark location into result
    *dest = location;
    ++dest;
    }
}
