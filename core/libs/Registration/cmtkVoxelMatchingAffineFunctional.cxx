/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkVoxelMatchingAffineFunctional.h>

#ifdef CMTK_BUILD_SMP
#  include <cmtkThreads.h>
#endif

#include <cmtkVoxelMatchingAffineFunctionalTemplate.h>
#include <cmtkVoxelMatchingMutInf.h>
#include <cmtkVoxelMatchingNormMutInf.h>
#include <cmtkVoxelMatchingCorrRatio.h>
#include <cmtkVoxelMatchingMeanSquaredDifference.h>
#include <cmtkVoxelMatchingCrossCorrelation.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

VoxelMatchingAffineFunctional* 
VoxelMatchingAffineFunctional
::Create
( const int metric, UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume, AffineXform::SmartPtr& affineXform )
{
  switch ( fltVolume->GetData()->GetDataClass() ) 
    {
    case DATACLASS_UNKNOWN :
    case DATACLASS_BINARY :
      // not a lot we can do here.
      return NULL;
    case DATACLASS_GREY :
      switch ( metric ) 
	{
	case 0:
	  return new VoxelMatchingAffineFunctionalTemplate< VoxelMatchingNormMutInf_Trilinear >( refVolume, fltVolume, affineXform );
	case 1:
	  return new VoxelMatchingAffineFunctionalTemplate<VoxelMatchingMutInf_Trilinear>( refVolume, fltVolume, affineXform );
	case 2:
	  return new VoxelMatchingAffineFunctionalTemplate<VoxelMatchingCorrRatio_Trilinear>( refVolume, fltVolume, affineXform );
	case 3:
	  return NULL; // masked nmi retired
	case 4:
	  return new VoxelMatchingAffineFunctionalTemplate<VoxelMatchingMeanSquaredDifference>( refVolume, fltVolume, affineXform );
	case 5:
	  return new VoxelMatchingAffineFunctionalTemplate<VoxelMatchingCrossCorrelation>( refVolume, fltVolume, affineXform );
	default:
	  break;
	}
      break;
    case DATACLASS_LABEL :
      switch ( metric ) 
	{
	case 0:
	  return new VoxelMatchingAffineFunctionalTemplate<VoxelMatchingNormMutInf_NearestNeighbor>( refVolume, fltVolume, affineXform );
	case 1:
	  return new VoxelMatchingAffineFunctionalTemplate<VoxelMatchingMutInf_NearestNeighbor>( refVolume, fltVolume, affineXform );
	case 2:
	  return new VoxelMatchingAffineFunctionalTemplate<VoxelMatchingCorrRatio_NearestNeighbor>( refVolume, fltVolume, affineXform );
	case 3:
	  return NULL; // masked nmi retired
	case 4:
	  return new VoxelMatchingAffineFunctionalTemplate<VoxelMatchingMeanSquaredDifference>( refVolume, fltVolume, affineXform );
	case 5:
	  return new VoxelMatchingAffineFunctionalTemplate<VoxelMatchingCrossCorrelation>( refVolume, fltVolume, affineXform );
	default:
	  break;
	}
      break;
    }
  
  return NULL;
}

} // namespace cmtk
