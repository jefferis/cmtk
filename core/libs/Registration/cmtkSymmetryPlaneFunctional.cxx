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

#include <cmtkSymmetryPlaneFunctional.h>

#include <cmtkVolumeAxesHash.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

SymmetryPlaneFunctional::SymmetryPlaneFunctional
( UniformVolume::SmartPtr& volume ) : m_Volume( NULL )
{
  this->SetVolume( volume );
  
  m_Metric = new VoxelMatchingNormMutInf<>( volume, volume );
}

SymmetryPlaneFunctional::SymmetryPlaneFunctional
( UniformVolume::SmartPtr& volume, 
  const Types::DataItemRange& valueRange )
  : m_Volume( NULL )
{
  this->SetVolume( volume );
  
  m_Metric = new VoxelMatchingNormMutInf<>( volume, volume, valueRange, valueRange );
}

Types::Coordinate 
SymmetryPlaneFunctional::GetParamStep 
( const size_t idx, const Types::Coordinate mmStep ) 
  const
{
  switch ( idx ) 
    {
    // plane offset is a translation
    case 0:
      return mmStep;
      // the other two parameters are rotations
    case 1:
    case 2:
      return mmStep / sqrt( MathUtil::Square( 0.5 * m_Volume->Size[0] ) + MathUtil::Square( 0.5 * m_Volume->Size[1] ) + MathUtil::Square( 0.5 * m_Volume->Size[2] ) ) * 90/M_PI;
    }
  return mmStep;
}

SymmetryPlaneFunctional::ReturnType
SymmetryPlaneFunctional::Evaluate()
{
  const VolumeAxesHash gridHash( *m_Volume, this->m_ParametricPlane, m_Volume->GetDelta() );
  const Vector3D *HashX = gridHash[0], *HashY = gridHash[1], *HashZ = gridHash[2];

  Vector3D pFloating;
    
  m_Metric->Reset();
    
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
	
	// Tell us whether the current location is still within the model
	// volume and get the respective voxel.
	if ( m_Volume->FindVoxelByIndex( pFloating, fltIdx, fltFrac ) ) 
	  {
	  // Compute data index of the model voxel in the model volume.
	  int offset = fltIdx[0] + DimsX * (fltIdx[1] + DimsY * fltIdx[2]);
	  
	  // Continue metric computation.
	  m_Metric->Proceed( (int) r, offset, fltFrac );
	  }
	}
      }
    }
  
  return m_Metric->Get();
}

} // namespace cmtk
