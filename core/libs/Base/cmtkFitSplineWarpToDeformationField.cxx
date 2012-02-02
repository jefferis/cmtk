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

#include "cmtkFitSplineWarpToDeformationField.h"

#include <Base/cmtkRegionIndexIterator.h>

#include <System/cmtkDebugOutput.h>

cmtk::FitSplineWarpToDeformationField::FitSplineWarpToDeformationField( DeformationField::SmartConstPtr dfield, const bool absolute ) 
  : m_DFieldIsAbsolute( absolute ), 
    m_DeformationField( dfield ),
    m_DeformationFieldFOV( dfield->m_Offset, dfield->m_Domain )
{  
}

cmtk::DataGrid::RegionType
cmtk::FitSplineWarpToDeformationField::GetDeformationGridRange( const UniformVolume::CoordinateRegionType& region ) const
{
  cmtk::DataGrid::IndexType regionFrom, regionTo;

  regionFrom = ComponentDivide( region.From() - this->m_DeformationField->m_Offset, this->m_DeformationField->m_Spacing );
  regionTo = ComponentDivide( region.To() - this->m_DeformationField->m_Offset, this->m_DeformationField->m_Spacing );

  regionFrom.AddScalar( 1 ); // to compensate for float-to-int truncation
  /// regionTo.AddScalar( 1 ); // necessary to convert to for() range, but COMMENT OUT to compensate for float-to-int truncation
  
  return cmtk::DataGrid::RegionType( regionFrom, regionTo );
}

void
cmtk::FitSplineWarpToDeformationField::ComputeResiduals( const SplineWarpXform& splineWarp )
{
  const DataGrid::IndexType dims = this->m_DeformationField->m_Dims;

  this->m_Residuals.resize( dims.Product() );

  size_t ofs = 0;
  for ( int z = 0; z < dims[2]; ++z )
    {
    for ( int y = 0; y < dims[1]; ++y )
      {
      for ( int x = 0; x < dims[0]; ++x, ++ofs )
	{
	this->m_Residuals[ofs] = this->m_DeformationField->GetTransformedGrid( x, y, z ) - splineWarp.GetTransformedGrid( x, y, z );
	if ( !this->m_DFieldIsAbsolute )
	  this->m_Residuals[ofs] += this->m_DeformationField->GetOriginalControlPointPosition( x, y, z );
	}
      }
    }
}

cmtk::SplineWarpXform::SmartPtr 
cmtk::FitSplineWarpToDeformationField::Fit( const SplineWarpXform::ControlPointIndexType& finalDims, const int nLevels )
{
  // compute initial dims based on nInitial = (nFinal-1) * 2^(nLevels-1) + 1
  SplineWarpXform::ControlPointIndexType initialDims = finalDims;
  initialDims.AddScalar( -1 );
  initialDims *= (1 << (nLevels-1));
  initialDims.AddScalar( +1 );

  // initialize B-spline transformation
  SplineWarpXform* splineWarp = new SplineWarpXform( this->m_DeformationField->m_Domain, initialDims, CoordinateVector::SmartPtr::Null(), AffineXform::SmartPtr( new AffineXform ) );
  
  this->FitSpline( *splineWarp, nLevels );
  
  return cmtk::SplineWarpXform::SmartPtr( splineWarp );
}

cmtk::SplineWarpXform::SmartPtr 
cmtk::FitSplineWarpToDeformationField::Fit( const Types::Coordinate finalSpacing, const int nLevels )
{
  // compute the start spacing of multi-level approximation by doubling final spacing until user-defined initial spacing is exceeded.
  Types::Coordinate spacing = finalSpacing * (1 << (nLevels-1));

  // initialize B-spline transformation
  SplineWarpXform* splineWarp = new SplineWarpXform( this->m_DeformationField->m_Domain, spacing, AffineXform::SmartPtr( new AffineXform ) );

  this->FitSpline( *splineWarp, nLevels );
  
  return cmtk::SplineWarpXform::SmartPtr( splineWarp );
}

void 
cmtk::FitSplineWarpToDeformationField::FitSpline( SplineWarpXform& splineWarp, const int nLevels )
{
  // loop until final control point spacing
  for ( int level = 0; level < nLevels; ++level )
    {
    DebugOutput( 5 ) << "Multi-resolution spline fitting level " << level+1 << " out of " << nLevels << "\n";

    // refine control point grid unless this is first iteration
    if ( level )
      {
      splineWarp.Refine();
      }

    // compute residuals
    splineWarp.RegisterVolumePoints( this->m_DeformationField->m_Dims, this->m_DeformationField->m_Spacing, this->m_DeformationField->m_Offset );

    this->ComputeResiduals( splineWarp );

    // loop over all control points to compute deltas as the spline coefficients that fit current residuals
    std::vector< FixedVector<3,Types::Coordinate> > delta( splineWarp.m_NumberOfControlPoints );

    const WarpXform::ControlPointRegionType cpRegion = splineWarp.GetAllControlPointsRegion();
    size_t cp = 0;
    for ( RegionIndexIterator<WarpXform::ControlPointRegionType> cpIt( cpRegion); cpIt != cpIt.end(); ++cpIt, ++cp )
      {
      // volume of influence for the current control point
      const DataGrid::RegionType voi = this->GetDeformationGridRange( splineWarp.GetVolumeOfInfluence( 3 * cp, this->m_DeformationFieldFOV, 0 /*fastMode=off*/ ) );
      
      // iterate over all voxels influenced by current control point.
      Types::Coordinate normalize = 0;
      for ( RegionIndexIterator<DataGrid::RegionType> it( voi ); it != it.end(); ++it )
	{
	const DataGrid::IndexType idx = it.Index();

	// Enumerator of Eq. (8) - this is a vector
	FixedVector<3,Types::Coordinate> ePc = this->m_Residuals[this->m_DeformationField->GetOffsetFromIndex( idx )/3]; // S_c(u1...un)
	Types::Coordinate pSquares = 1; // this is the product over the B^2-s in Eq. (11)
	for ( int axis = 0; axis < 3; ++axis )
	  {
	  // relative index of spline function for current pixel relative to current control point
	  const int relIdx = cpIt.Index()[axis] - splineWarp.m_GridIndexes[axis][idx[axis]];
	  // sanity range checks
	  assert( (relIdx >= 0) && (relIdx < 4) );

	  ePc *= splineWarp.m_GridSpline[axis][4*it.Index()[axis]+relIdx];
	  pSquares *= MathUtil::Square( splineWarp.m_GridSpline[axis][4*it.Index()[axis]+relIdx] );
	  }

	// Denominator of Eq. (8) - this is a scalar
	Types::Coordinate dPc = 0;

	const DataGrid::RegionType neighborhood( DataGrid::IndexType::Init( 0 ), DataGrid::IndexType::Init( 4 ) );
	for ( RegionIndexIterator<DataGrid::RegionType> nIt( neighborhood ); nIt != nIt.end(); ++nIt )
	  {
	  Types::Coordinate prod = 1;
	  for ( int axis = 0; axis < 3; ++axis )
	    {
	    prod *= MathUtil::Square( splineWarp.m_GridSpline[axis][it.Index()[axis]+nIt.Index()[axis]] );
	    }
	  
	  dPc += prod;
	  }

	// Eq. (11)
	delta[cp] += (pSquares / dPc) * ePc;

	normalize += pSquares;
	}
      
      // Eq. (11) denominator
      delta[cp] /= normalize;
      }
    
    // apply delta
    for ( size_t cp = 0; cp < splineWarp.m_NumberOfControlPoints; ++cp )
      {
      splineWarp.SetShiftedControlPointPositionByOffset( splineWarp.GetShiftedControlPointPositionByOffset( cp ) + delta[cp], cp );
      }
    }
}
