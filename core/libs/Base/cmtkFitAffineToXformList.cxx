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

#include "cmtkFitAffineToXformList.h"

#include <Base/cmtkRegionIndexIterator.h>

#include <System/cmtkDebugOutput.h>

cmtk::AffineXform::SmartPtr 
cmtk::FitAffineToXformList::Fit()
{
  // first, get the centroids in "from" and "to" space
  cmtk::FixedVector<3,cmtk::Types::Coordinate> cFrom( 0.0 );
  cmtk::FixedVector<3,cmtk::Types::Coordinate> cTo( 0.0 );

  size_t numberOfValidPixels = 0;
  
  size_t ofs = 0;
  for ( RegionIndexIterator<DataGrid::RegionType> it = this->m_XformField.GetWholeImageRegion(); it != it.end(); ++it, ++ofs )
    {
    if ( this->m_XformValidAt[ofs] )
      {
      cFrom += this->m_XformField.GetGridLocation( it.Index() );
      cTo += this->m_XformField[ofs];
      ++numberOfValidPixels;
      }
    }
  
  cFrom /= numberOfValidPixels;
  cTo /= numberOfValidPixels;
  
  // now get the transformation matrix for rotation, scale, shear, using the previously computed centroids for reference
  Matrix3x3<Types::Coordinate> matrix = this->GetMatrix( cFrom, cTo );
  
  // put everything together
  AffineXform::MatrixType matrix4x4( matrix );
  AffineXform::SmartPtr result( new AffineXform( matrix4x4 ) );
  result->SetTranslation( (cTo - cFrom) );
  result->SetCenter( cFrom );
  
  return result;
}

cmtk::Matrix3x3<cmtk::Types::Coordinate> 
cmtk::FitAffineToXformList::GetMatrix( const cmtk::FixedVector<3,cmtk::Types::Coordinate>& cFrom, const cmtk::FixedVector<3,cmtk::Types::Coordinate>& cTo )
{
  Matrix3x3<Types::Coordinate> txT = Matrix3x3<Types::Coordinate>::Zero(); // "t" is the 3xN matrix of transformation vectors (after removing global translation) at the control points
  Matrix3x3<Types::Coordinate> xxT = Matrix3x3<Types::Coordinate>::Zero(); // "x" is the 3xN matrix of control point grid coordinates
  
  // build the two 3x3 matrices of (t*xT)(x*xT) on the fly.
  size_t ofs = 0;
  for ( RegionIndexIterator<DataGrid::RegionType> it = this->m_XformField.GetWholeImageRegion(); it != it.end(); ++it, ++ofs )
    {
    if ( this->m_XformValidAt[ofs] )
      {
      const cmtk::FixedVector<3,cmtk::Types::Coordinate> x = this->m_XformField.GetGridLocation( it.Index() ) - cFrom;
      const cmtk::FixedVector<3,cmtk::Types::Coordinate> t = this->m_XformField[ofs] - cTo;
      
      for ( size_t j = 0; j < 3; ++j )
	{
	for ( size_t i = 0; i < 3; ++i )
	  {
	  txT[i][j] += t[j] * x[i];
	  xxT[i][j] += x[j] * x[i];
	  }
	}
      }
    }
  
  return (xxT.GetInverse()*txT);
}
