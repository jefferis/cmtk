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

#include "cmtkFitAffineToWarpXform.h"

#include <Base/cmtkRegionIndexIterator.h>

#include <System/cmtkDebugOutput.h>

cmtk::FitAffineToWarpXform::FitAffineToWarpXform( WarpXform::SmartConstPtr warp ) 
  : m_WarpXform( warp )
{  
}

cmtk::AffineXform::SmartPtr 
cmtk::FitAffineToWarpXform::Fit()
{
  const cmtk::FixedVector<3,cmtk::Types::Coordinate> xlate = Self::GetMeanTranslation( *(this->m_WarpXform) );

  Matrix3x3<Types::Coordinate> matrix = Self::GetMatrix( *(this->m_WarpXform), xlate );
  
  AffineXform::MatrixType matrix4x4( matrix );
  AffineXform::SmartPtr result( new AffineXform( matrix4x4 ) );
  result->SetTranslation( xlate );
  
  return result;
}

cmtk::FixedVector<3,cmtk::Types::Coordinate> 
cmtk::FitAffineToWarpXform::GetMeanTranslation( const WarpXform& warpXform )
{
  cmtk::FixedVector<3,cmtk::Types::Coordinate> delta( cmtk::FixedVector<3,cmtk::Types::Coordinate>::Init( 0 ) );
  
  size_t ofs = 0;
  for ( RegionIndexIterator<WarpXform::ControlPointRegionType> it = warpXform.GetAllControlPointsRegion(); it != it.end(); ++it, ++ofs )
    {
    delta += warpXform.GetDeformedControlPointPosition( it.Index()[0], it.Index()[1], it.Index()[2] ) - warpXform.GetOriginalControlPointPositionByOffset( ofs );
    }

  return (delta /= warpXform.GetAllControlPointsRegion().Size());
}

cmtk::Matrix3x3<cmtk::Types::Coordinate> 
cmtk::FitAffineToWarpXform::GetMatrix( const WarpXform& warpXform, const cmtk::FixedVector<3,cmtk::Types::Coordinate>& xlate )
{
  Matrix3x3<Types::Coordinate> txT; // "t" is the 3xN matrix of transformation vectors (after removing global translation) at the control points
  Matrix3x3<Types::Coordinate> xxT; // "x" is the 3xN matrix of control point grid coordinates

  txT.Fill( 0.0 );
  xxT.Fill( 0.0 );

  size_t ofs = 0;
  for ( RegionIndexIterator<WarpXform::ControlPointRegionType> it = warpXform.GetAllControlPointsRegion(); it != it.end(); ++it, ++ofs )
    {
    const cmtk::FixedVector<3,cmtk::Types::Coordinate> x = warpXform.GetOriginalControlPointPositionByOffset( ofs );
    const cmtk::FixedVector<3,cmtk::Types::Coordinate> t = warpXform.GetDeformedControlPointPosition( it.Index()[0], it.Index()[1], it.Index()[2] ) - x - xlate;

    for ( size_t j = 0; j < 3; ++j )
      {
      for ( size_t i = 0; i < 3; ++i )
	{
	txT[i][j] += t[i] * x[j];
	xxT[i][j] += x[i] * x[j];
	}
      }
    }  
  
  return (txT * xxT.Invert3x3());
}
