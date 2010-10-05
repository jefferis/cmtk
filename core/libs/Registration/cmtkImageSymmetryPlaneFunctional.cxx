/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkImageSymmetryPlaneFunctional.h"

#include <Base/cmtkTransformedVolumeAxes.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImageSymmetryPlaneFunctional::ImageSymmetryPlaneFunctional
( UniformVolume::SmartConstPtr& volume ) 
  : ImageSymmetryPlaneFunctionalBase( volume ),
    m_Metric( new ImagePairSimilarityMeasureMSD( this->m_Volume, this->m_Volume ) )
{
}

ImageSymmetryPlaneFunctional::ImageSymmetryPlaneFunctional
( UniformVolume::SmartConstPtr& volume, 
  const Types::DataItemRange& valueRange )
  : ImageSymmetryPlaneFunctionalBase( volume, valueRange ),
    m_Metric( new ImagePairSimilarityMeasureMSD( this->m_Volume, this->m_Volume ) )
{
}

ImageSymmetryPlaneFunctional::ReturnType
ImageSymmetryPlaneFunctional::Evaluate()
{
  const TransformedVolumeAxes gridHash( *m_Volume, this->m_ParametricPlane, m_Volume->Deltas().begin() );
  const Vector3D *HashX = gridHash[0], *HashY = gridHash[1], *HashZ = gridHash[2];

  Vector3D pFloating;
    
  Self::MetricType& metric = *m_Metric;
  metric.Reset();
    
  const DataGrid::IndexType& Dims = m_Volume->GetDims();
  const int DimsX = Dims[0], DimsY = Dims[1], DimsZ = Dims[2];

  int fltIdx[3];
  Types::Coordinate fltFrac[3];

  Vector3D planeStart, rowStart;

  int r = 0;
  for ( int pZ = 0; pZ<DimsZ; ++pZ ) 
    {
    planeStart = HashZ[pZ];
    
    for ( int pY = 0; pY<DimsY; ++pY ) 
      {
      (rowStart = planeStart) += HashY[pY];
      
      for ( int pX = 0; pX<DimsX; ++pX, ++r ) 
	{
	(pFloating = rowStart) += HashX[pX];
	
	// Is the current location still within the floating image, then get the respective voxel.
	if ( m_Volume->FindVoxelByIndex( pFloating, fltIdx, fltFrac ) )
	  {
	  // Continue metric computation.
	  metric.Increment( metric.GetSampleX( r ), metric.GetSampleY( fltIdx, fltFrac ) );
	  }
	}
      }
    }
  
  return metric.Get();
}

} // namespace cmtk
