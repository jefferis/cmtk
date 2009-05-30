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

#include <cmtkSplineWarpCongealingFunctional.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

void
SplineWarpCongealingFunctional::StaticThreadStorage::Initialize( const Self* This )
{
  const size_t numberOfXforms = This->m_XformVector.size();
  
  this->m_FPlus.resize( 3 * numberOfXforms );
  this->m_FMinus.resize( 3 * numberOfXforms );
  this->m_CountByParameterPlus.resize( 3 * numberOfXforms );
  this->m_CountByParameterMinus.resize( 3 * numberOfXforms );
  
  this->m_Xforms.resize( numberOfXforms );
  for ( size_t xi = 0; xi < numberOfXforms; ++xi )
    {
    this->m_Xforms[xi] = SplineWarpXform::SmartPtr( This->GetXformByIndex(xi)->Clone() );
    }
  
  this->m_VectorList.resize( This->m_MaximumNumberOfPixelsPerLineVOI );
  this->m_Count.resize( This->m_MaximumNumberOfPixelsPerLineVOI );
  this->m_Histogram.resize( This->m_MaximumNumberOfPixelsPerLineVOI );
  for ( size_t x = 0; x < This->m_MaximumNumberOfPixelsPerLineVOI; ++x )
    {
    this->m_Histogram[x].Resize( This->m_HistogramBins + 2 * This->m_HistogramKernelRadiusMax, false /*reset*/ );
    }    
}

} // namespace cmtk
