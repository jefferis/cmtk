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

#include "cmtkEchoPlanarUnwarpFunctional.h"

#include <Base/cmtkDataGrid.h>
#include <Base/cmtkSincInterpolator.h>
#include <Base/cmtkRegionIndexIterator.h>

#include <algorithm>

const int cmtk::EchoPlanarUnwarpFunctional::InterpolationKernelRadius = 3; 

void
cmtk::EchoPlanarUnwarpFunctional::ComputeDeformedImage( const UniformVolume& sourceImage, std::vector<Types::DataItem>& targetImageData, int direction )
{
  const Types::Coordinate pixelSize = sourceImage.Deltas()[this->m_PhaseEncodeDirection];

  const DataGrid::RegionType wholeImageRegion = sourceImage.CropRegion();
  for ( RegionIndexIterator<DataGrid::RegionType> it( wholeImageRegion ); it != it.end(); ++it )
    {
    DataGrid::IndexType idx = it.Index();
    const size_t i = sourceImage.GetOffsetFromIndex( idx );

    // first, get Jacobian for grid position
    targetImageData[i] = 1 + direction * this->GetPartialJacobian( idx );

    // now compute deformed position for interpolation
    const Types::Coordinate shift = direction*this->m_Deformation[i] / pixelSize;
    const Types::Coordinate position = shift + idx[this->m_PhaseEncodeDirection];
    
    idx[this->m_PhaseEncodeDirection] = static_cast<int>( floor( position ) );

    // multiple interpolated data onto previously set Jacobian
    targetImageData[i] *= this->Interpolate1D( sourceImage, idx, position - idx[this->m_PhaseEncodeDirection] );    
    }
}

cmtk::Types::DataItem 
cmtk::EchoPlanarUnwarpFunctional::Interpolate1D( const UniformVolume& sourceImage, const FixedVector<3,int>& baseIdx, const Types::Coordinate relative ) const
{
  FixedVector<3,int> idx = baseIdx;

  const int iFrom = -std::min( Self::InterpolationKernelRadius, idx[this->m_PhaseEncodeDirection] );
  const int iTo = std::max( Self::InterpolationKernelRadius, sourceImage.m_Dims[this->m_PhaseEncodeDirection] - idx[this->m_PhaseEncodeDirection] - 1 );

  idx[this->m_PhaseEncodeDirection] += iFrom;

  Types::DataItem value = 0;
  Types::Coordinate total = 0;

  for ( int i = iFrom; i < iTo; ++i, ++idx[this->m_PhaseEncodeDirection] )
    {
    const Types::Coordinate weight = Interpolators::CosineSinc<Self::InterpolationKernelRadius>::GetWeight( i, relative );
    value += weight * sourceImage.GetDataAt( sourceImage.GetOffsetFromIndex( idx ) );
    total += weight;
    }
  
  if ( total > 0 )
    return static_cast<Types::DataItem>( value / total );
  else
    return 0;
}

cmtk::Types::Coordinate 
cmtk::EchoPlanarUnwarpFunctional::GetPartialJacobian( const FixedVector<3,int>& baseIdx ) const
{
  cmtk::Types::Coordinate diff = 0;
  int normalize = 0;

  size_t offset = this->m_ImageFwd->GetOffsetFromIndex( baseIdx );
  if ( baseIdx[this->m_PhaseEncodeDirection] > 0 )
    {
    diff -= this->m_Deformation[ offset - this->m_ImageGrid->m_GridIncrements[this->m_PhaseEncodeDirection] ];
    ++normalize;
    }
  else
    {
    diff -= this->m_Deformation[offset];
    }

  if ( baseIdx[this->m_PhaseEncodeDirection] < this->m_ImageGrid->m_Dims[this->m_PhaseEncodeDirection]-1 )
    {
    diff += this->m_Deformation[ offset + this->m_ImageGrid->m_GridIncrements[this->m_PhaseEncodeDirection] ];
    ++normalize;
    }
  else
    {
    diff += this->m_Deformation[offset];
    }
  
  return diff / normalize;
}

void
cmtk::EchoPlanarUnwarpFunctional
::FunctionAndGradient
::Evaluate( const ap::real_1d_array& x, ap::real_value_type& f, ap::real_1d_array& g )
{
  this->m_Function->ComputeDeformedImage( *(this->m_Function->m_ImageFwd), this->m_Function->m_UnwarpImageFwd, +1 );
  this->m_Function->ComputeDeformedImage( *(this->m_Function->m_ImageRev), this->m_Function->m_UnwarpImageRev, -1 );

  const size_t nPixels = this->m_Function->m_ImageGrid->GetNumberOfPixels();
  ap::real_value_type msd = 0;
  for ( size_t px = 0; px < nPixels; ++px )
    {
    msd += MathUtil::Square( this->m_Function->m_UnwarpImageFwd[px] - this->m_Function->m_UnwarpImageRev[px] );
    }

  f = msd / nPixels;

  for ( size_t px = 0; px < nPixels; ++px )
    {
    
    }
}
