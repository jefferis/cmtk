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

#include "cmtkLeastSquaresPolynomialIntensityBiasField.h"

#include <Base/cmtkPolynomial.h>
#include <Base/cmtkRegionIndexIterator.h>
#include <Base/cmtkMatrix.h>
#include <Base/cmtkLeastSquares.h>

namespace
cmtk
{

LeastSquaresPolynomialIntensityBiasField::LeastSquaresPolynomialIntensityBiasField( const UniformVolume& image, const std::vector<bool>& mask, const int degree )
{
  const UniformVolume::CoordinateVectorType center = image.GetCenterCropRegion();

  // first, compute average intensity over masked region
  Types::DataItem avg = 0;
  size_t nPixelsMask = 0;

  const DataGrid::RegionType region = image.GetWholeImageRegion();
  for ( RegionIndexIterator<DataGrid::RegionType> it( region ); it != it.end(); ++it )
    {
    const size_t ofs = image.GetOffsetFromIndex( it.Index() );
    if ( mask[ofs] )
      {
      avg += fabs( image.GetDataAt( ofs ) );
      ++nPixelsMask;
      }
    }

  if ( !nPixelsMask )
    throw Self::EmptyMaskException();

  avg /= nPixelsMask;

  // set up least-squares problem
  size_t nVars = 0;
  switch ( degree )
    {
    case 1: nVars = Polynomial<1,Types::DataItem>::NumberOfMonomials; break;
    case 2: nVars = Polynomial<2,Types::DataItem>::NumberOfMonomials; break;
    case 3: nVars = Polynomial<3,Types::DataItem>::NumberOfMonomials; break;
    case 4: nVars = Polynomial<4,Types::DataItem>::NumberOfMonomials; break;
    }

  std::vector<Types::DataItem> dataVector( nPixelsMask );
  Matrix2D<Types::DataItem> uMatrix( nPixelsMask, nVars );

  size_t cntPx = 0;
  for ( RegionIndexIterator<DataGrid::RegionType> it( region ); it != it.end(); ++it )
    {
    const size_t ofs = image.GetOffsetFromIndex( it.Index() );

    if ( mask[ofs] )
      {
      const UniformVolume::CoordinateVectorType xyz = ComponentDivide( image.GetGridLocation( it.Index() ) - center, image.Size );    
      dataVector[cntPx] = image.GetDataAt( ofs ) / avg - 1.0;
      for ( size_t n = 0; n < nVars; ++n )
	{
	uMatrix[cntPx][n] = Polynomial<4,Types::DataItem>::EvaluateMonomialAt( n, xyz[0], xyz[1], xyz[2] );
	}
      ++cntPx;
      }
    }

  // solve least-squares problem
  const std::vector<Types::DataItem> params = LeastSquares<Types::DataItem>( uMatrix ).Solve( dataVector );

  // apply solution
  this->m_CorrectedData = TypedArray::Create( image.GetData()->GetType(), image.GetNumberOfPixels() );
  this->m_BiasData = TypedArray::Create( TYPE_ITEM, image.GetNumberOfPixels() );

  for ( RegionIndexIterator<DataGrid::RegionType> it( region ); it != it.end(); ++it )
    {
    const size_t ofs = image.GetOffsetFromIndex( it.Index() );

    const UniformVolume::CoordinateVectorType xyz = ComponentDivide( image.GetGridLocation( it.Index() ) - center, image.Size );    
    
    Types::DataItem bias = 1.0;
    for ( size_t n = 0; n < nVars; ++n )
      {
      bias += params[n] * Polynomial<4,Types::DataItem>::EvaluateMonomialAt( n, xyz[0], xyz[1], xyz[2] );
      }
    
    this->m_BiasData->Set( bias, ofs );

    Types::DataItem value;
    if ( image.GetData()->Get( value, ofs ) )
      {
      this->m_CorrectedData->Set( value / bias, ofs );
      }
    else
      {
      this->m_CorrectedData->SetPaddingAt( ofs );
      }
    }
}

} // namespace cmtk
