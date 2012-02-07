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

#pragma omp parallel for
  for ( int z = 0; z < dims[2]; ++z )
    {
    size_t ofs = z * dims[0] * dims[1];

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
  // we may have to adjust nLevels downwards
  int numberOfLevels = nLevels;
  
  SplineWarpXform::ControlPointIndexType initialDims = finalDims;
  for ( int level = 1; level < numberOfLevels; ++level )
    {
    if ( (initialDims[0]&1) && (initialDims[1]&1) && (initialDims[2]&1) && // check that all dims are still odd numbers
	 (initialDims.MinValue()>4) ) // check that at least 4 control points will be left in each dimension
      {
      // apply the inverse of the refinement formula used in SplineWarpXform::Refine()
      initialDims.AddScalar( +3 );
      initialDims /= 2;
      }
    else
      {
      numberOfLevels = level;

      DebugOutput(2) << "INFO: adjusted number of levels to " << numberOfLevels << " from " << nLevels << " to ensure sufficient number of control points\n";
      }
    }

  // initialize B-spline transformation
  SplineWarpXform* splineWarp = new SplineWarpXform( this->m_DeformationField->m_Domain, initialDims, CoordinateVector::SmartPtr::Null(), AffineXform::SmartPtr( new AffineXform ) );
  
  this->FitSpline( *splineWarp, numberOfLevels );
  
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

    DebugOutput( 6 ) << "  Control point grid is " << splineWarp.m_Dims[0] << "x" << splineWarp.m_Dims[1] << "x" << splineWarp.m_Dims[2] << "\n";

    // compute residuals
    splineWarp.RegisterVolumePoints( this->m_DeformationField->m_Dims, this->m_DeformationField->m_Spacing, this->m_DeformationField->m_Offset );

    this->ComputeResiduals( splineWarp );

    // loop over all control points to compute deltas as the spline coefficients that fit current residuals
    std::vector< FixedVector<3,Types::Coordinate> > delta( splineWarp.m_NumberOfControlPoints, FixedVector<3,Types::Coordinate>( FixedVector<3,Types::Coordinate>::Init( 0.0 ) ) );

    const WarpXform::ControlPointRegionType cpRegionAll = splineWarp.GetAllControlPointsRegion();
#ifndef _OPENMP
    const WarpXform::ControlPointRegionType cpRegion = cpRegionAll;
#else // _OPENMP
    const size_t maxIdx = (cpRegionAll.To()-cpRegionAll.From()).MaxIndex();
    const int sliceFrom = cpRegionAll.From()[maxIdx];
    const int sliceTo = cpRegionAll.To()[maxIdx];
#pragma omp parallel for
    for ( int slice = sliceFrom; slice < sliceTo; ++slice )
      {
      WarpXform::ControlPointRegionType cpRegion = cpRegionAll;
      cpRegion.From()[maxIdx] = slice;
      cpRegion.To()[maxIdx] = slice+1;
#endif
      for ( RegionIndexIterator<WarpXform::ControlPointRegionType> cpIt( cpRegion); cpIt != cpIt.end(); ++cpIt )
	{
	const size_t cp = splineWarp.GetOffsetFromIndex( cpIt.Index() ) / 3;
	
	// volume of influence for the current control point
	const DataGrid::RegionType voi = this->GetDeformationGridRange( splineWarp.GetVolumeOfInfluence( 3 * cp, this->m_DeformationFieldFOV, false /*fastMode=off*/ ) );
	
	// iterate over all voxels influenced by current control point.
	Types::Coordinate normalize = 0;
	for ( RegionIndexIterator<DataGrid::RegionType> it( voi ); it != it.end(); ++it )
	  {
	  const DataGrid::IndexType idx = it.Index();
	  
	  // Enumerator of Eq. (8) - this is a vector
	  Types::Coordinate pB = 1; // this is the product over the B in enumerator of Eq. (8)
	  for ( int axis = 0; axis < 3; ++axis )
	    {
	    // relative index of spline function for current pixel relative to current control point
	    const int relIdx = cpIt.Index()[axis] - splineWarp.m_GridIndexes[axis][idx[axis]];
	    // sanity range checks
	    assert( (relIdx >= 0) && (relIdx < 4) );
	    
	    pB *= splineWarp.m_GridSpline[axis][4*it.Index()[axis]+relIdx];
	    }
	  
	  // Denominator of Eq. (8) - this is a scalar
	  Types::Coordinate dPc = 0;
	  
	  const DataGrid::RegionType neighborhood( DataGrid::IndexType::Init( 0 ), DataGrid::IndexType::Init( 4 ) );
	  for ( RegionIndexIterator<DataGrid::RegionType> nIt( neighborhood ); nIt != nIt.end(); ++nIt )
	    {
	    Types::Coordinate prod = 1;
	    for ( int axis = 0; axis < 3; ++axis )
	      {
	      prod *= splineWarp.m_GridSpline[axis][(4*it.Index()[axis])+nIt.Index()[axis]];
	      }
	    
	    dPc += MathUtil::Square( prod );
	    }
	  
	  // Eq. (11)
	  const Types::Coordinate pB2 = MathUtil::Square( pB );
	  delta[cp] += pB2 * (pB / dPc) * this->m_Residuals[this->m_DeformationField->GetOffsetFromIndex( idx )/3]; // S_c(u1...un)
	  normalize += pB2;
	  }
	
	// Eq. (11) denominator
	delta[cp] /= normalize;
	}
#ifdef _OPENMP
      }
#endif
    
    // apply delta
    for ( size_t cp = 0; cp < splineWarp.m_NumberOfControlPoints; ++cp )
      {
      splineWarp.SetShiftedControlPointPositionByOffset( splineWarp.GetShiftedControlPointPositionByOffset( cp ) + delta[cp], cp );
      }
    }
}
