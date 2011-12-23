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
#include <Base/cmtkUniformVolumeFilter.h>

#include <System/cmtkDebugOutput.h>
#include <System/cmtkProgress.h>

#include <IO/cmtkVolumeIO.h>

#include <algorithm>

const int cmtk::EchoPlanarUnwarpFunctional::InterpolationKernelRadius = 3; 

cmtk::EchoPlanarUnwarpFunctional::EchoPlanarUnwarpFunctional
( UniformVolume::SmartConstPtr& imageFwd, UniformVolume::SmartConstPtr& imageRev, const byte phaseEncodeDirection )
  : m_ImageGrid( imageFwd->CloneGrid() ), 
    m_ImageFwd( imageFwd ), 
    m_ImageRev( imageRev ), 
    m_SmoothImageFwd( imageFwd ), 
    m_SmoothImageRev( imageRev ), 
    m_PhaseEncodeDirection( phaseEncodeDirection )
{
  this->m_Deformation.setbounds( 1, this->m_ImageGrid->GetNumberOfPixels() );
  for ( size_t i = 1; i < 1+this->m_ImageGrid->GetNumberOfPixels(); ++i )
    this->m_Deformation(i) = 0.0;

  this->m_UnwarpImageFwd.resize( this->m_ImageGrid->GetNumberOfPixels() );
  this->m_UnwarpImageRev.resize( this->m_ImageGrid->GetNumberOfPixels() );

  this->m_CorrectedImageFwd.resize( this->m_ImageGrid->GetNumberOfPixels() );
  this->m_CorrectedImageRev.resize( this->m_ImageGrid->GetNumberOfPixels() );
}

void
cmtk::EchoPlanarUnwarpFunctional::SetSmoothingKernelWidth( const Units::GaussianSigma& sigma, const Types::Coordinate maxError )
{
  if ( sigma.Value() > 0 )
    {
    {
    UniformVolumeFilter filterFwd( this->m_ImageFwd );
    UniformVolume::SmartPtr smooth = UniformVolume::SmartPtr( this->m_ImageGrid->CloneGrid() );
    smooth->SetData( filterFwd.GetDataGaussFiltered1D( this->m_PhaseEncodeDirection, sigma, maxError ) );
    this->m_SmoothImageFwd = smooth;

    VolumeIO::Write( *(this->m_SmoothImageFwd), "smoothF.nii" );
    }

    {
    UniformVolumeFilter filterRev( this->m_ImageRev );
    UniformVolume::SmartPtr smooth = UniformVolume::SmartPtr( this->m_ImageGrid->CloneGrid() );
    smooth->SetData( filterRev.GetDataGaussFiltered1D( this->m_PhaseEncodeDirection, sigma, maxError ) );
    this->m_SmoothImageRev = smooth;

    VolumeIO::Write( *(this->m_SmoothImageRev), "smoothR.nii" );
    }    
    }
  else
    {
    this->m_SmoothImageFwd = this->m_ImageFwd;
    this->m_SmoothImageRev = this->m_ImageRev;
    }
}


void
cmtk::EchoPlanarUnwarpFunctional::MakeGradientImage( const ap::real_1d_array& u, const int direction, const UniformVolume& sourceImage, std::vector<Types::DataItem>& gradientImageData )
{
  DebugOutput( 9 ) << "Making gradient image\n";

  gradientImageData.resize( sourceImage.GetNumberOfPixels() );

  const DataGrid::RegionType wholeImageRegion = sourceImage.GetWholeImageRegion();

#ifndef _OPENMP
  const DataGrid::RegionType region = wholeImageRegion;
#else // _OPENMP
#pragma omp parallel for
  for ( int slice = wholeImageRegion.From()[2]; slice < wholeImageRegion.To()[2]; ++slice )
    {
    DataGrid::RegionType region = wholeImageRegion;
    region.From()[2] = slice;
    region.To()[2] = slice+1;
#endif
    for ( RegionIndexIterator<DataGrid::RegionType> it( region ); it != it.end(); ++it )
      {
      DataGrid::IndexType idx = it.Index();
      const size_t i = sourceImage.GetOffsetFromIndex( idx );
      
      // apply deformation
      const Types::Coordinate shift = direction * u(1+i) + idx[this->m_PhaseEncodeDirection];

      Types::Coordinate position = shift + 0.5;      
      idx[this->m_PhaseEncodeDirection] = static_cast<int>( floor( position ) );
      gradientImageData[i] = this->Interpolate1D( sourceImage, idx, position -  idx[this->m_PhaseEncodeDirection] );
      
      position = shift - 0.5;      
      idx[this->m_PhaseEncodeDirection] = static_cast<int>( floor( position ) );
      gradientImageData[i] -= this->Interpolate1D( sourceImage, idx, position -  idx[this->m_PhaseEncodeDirection] );
      }
#ifdef _OPENMP
    }
#endif
}

cmtk::Types::DataItem 
cmtk::EchoPlanarUnwarpFunctional::Interpolate1D( const UniformVolume& sourceImage, const FixedVector<3,int>& baseIdx, const Types::Coordinate relative ) const
{
  FixedVector<3,int> idx = baseIdx;

  const int maxIdx = sourceImage.m_Dims[this->m_PhaseEncodeDirection] - 1;

  const int iFrom = -std::min( Self::InterpolationKernelRadius, idx[this->m_PhaseEncodeDirection] );
  const int iTo = std::min( Self::InterpolationKernelRadius, maxIdx - idx[this->m_PhaseEncodeDirection] );
  
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

void
cmtk::EchoPlanarUnwarpFunctional::ComputeDeformedImage( const ap::real_1d_array& u, int direction, const UniformVolume& sourceImage, std::vector<Types::DataItem>& targetUnwarpData, std::vector<Types::DataItem>& targetCorrectedData )
{
  DebugOutput( 9 ) << "Computing deformed image\n";

  const DataGrid::RegionType wholeImageRegion = sourceImage.GetWholeImageRegion();

#ifndef _OPENMP
  const DataGrid::RegionType region = wholeImageRegion;
#else // _OPENMP
#pragma omp parallel for
  for ( int slice = wholeImageRegion.From()[2]; slice < wholeImageRegion.To()[2]; ++slice )
    {
    DataGrid::RegionType region = wholeImageRegion;
    region.From()[2] = slice;
    region.To()[2] = slice+1;
#endif
    for ( RegionIndexIterator<DataGrid::RegionType> it( region ); it != it.end(); ++it )
      {
      DataGrid::IndexType idx = it.Index();
      const size_t i = sourceImage.GetOffsetFromIndex( idx );
      
      // first, get Jacobian for grid position
      targetCorrectedData[i] = 1 + direction * this->GetPartialJacobian( u, idx );
      
      // now compute deformed position for interpolation
      const Types::Coordinate shift = direction * u(1+i);
      const Types::Coordinate position = shift + idx[this->m_PhaseEncodeDirection];
      
      idx[this->m_PhaseEncodeDirection] = static_cast<int>( floor( position ) );
      
      // multiple interpolated data onto previously set Jacobian
      targetCorrectedData[i] *= (targetUnwarpData[i] = this->Interpolate1D( sourceImage, idx, position - idx[this->m_PhaseEncodeDirection] ));
      }
#ifdef _OPENMP
    }
#endif
}

cmtk::Types::Coordinate 
cmtk::EchoPlanarUnwarpFunctional::GetPartialJacobian( const ap::real_1d_array& u, const FixedVector<3,int>& baseIdx ) const
{
  cmtk::Types::Coordinate diff = 0;
  int normalize = 0;

  size_t offset = this->m_ImageGrid->GetOffsetFromIndex( baseIdx );
  if ( baseIdx[this->m_PhaseEncodeDirection] > 0 )
    {
    diff -= u( 1 + offset - this->m_ImageGrid->m_GridIncrements[this->m_PhaseEncodeDirection] );
    ++normalize;
    }
  else
    {
    diff -= u( 1 + offset );
    }

  if ( baseIdx[this->m_PhaseEncodeDirection] < this->m_ImageGrid->m_Dims[this->m_PhaseEncodeDirection]-1 )
    {
    diff += u( 1 + offset + this->m_ImageGrid->m_GridIncrements[this->m_PhaseEncodeDirection] );
    ++normalize;
    }
  else
    {
    diff += u( 1 + offset );
    }
  
  return diff / normalize;
}

void
cmtk::EchoPlanarUnwarpFunctional::Optimize( const int numberOfIterations, const Units::GaussianSigma& smoothMax, const Units::GaussianSigma& smoothDiff )
{
  const int numberOfPixels = this->m_ImageGrid->GetNumberOfPixels();
  
  ap::integer_1d_array nbd;
  nbd.setbounds( 1, numberOfPixels );
  for ( int i = 1; i <= numberOfPixels; ++i )
    {
    nbd(i) = 0;
    }

  // dummy array for unused constraints
  ap::real_1d_array dummy;
  
  for ( Units::GaussianSigma smoothness = smoothMax; smoothness.Value() >= 0; smoothness = smoothness - smoothDiff )
    {
    DebugOutput( 4 ) << "Setting image smoothing kernel sigma=" << smoothness.Value() << "\n";
    this->SetSmoothingKernelWidth( smoothness );
    
    Progress::Begin( 0, numberOfIterations, 1, "EPI Unwarping" );
    
    int info;
    Self::FunctionAndGradient functionAndGradient( this );
    ap::lbfgsbminimize( &(functionAndGradient), numberOfPixels, 5, this->m_Deformation, 1e-10 /*epsg*/, 1e-10 /*epsf*/, 1e-10 /*epsx*/, numberOfIterations, nbd, dummy, dummy, info );
    
    Progress::Done();
    
    if ( info < 0 )
      StdErr << "ERROR: lbfgsbminimize returned status code " << info << "\n";
    }

  // update corrected images using unsmoothed images with final deformation
  this->ComputeDeformedImage( this->m_Deformation, +1, *(this->m_ImageFwd), this->m_UnwarpImageFwd, this->m_CorrectedImageFwd );
  this->ComputeDeformedImage( this->m_Deformation, -1, *(this->m_ImageRev), this->m_UnwarpImageRev, this->m_CorrectedImageRev );
}


void
cmtk::EchoPlanarUnwarpFunctional
::FunctionAndGradient
::Evaluate( const ap::real_1d_array& x, ap::real_value_type& f, ap::real_1d_array& g )
{
  FunctionType& function = *(this->m_Function);

  // reset gradient vector
  for ( size_t i = g.getlowbound(); i <= g.gethighbound(); ++i )
    g(i) = 0.0;

  const UniformVolume& sourceImage = *(function.m_ImageGrid);
  const DataGrid::RegionType wholeImageRegion = sourceImage.GetWholeImageRegion();
  
  function.ComputeDeformedImage( x, +1, *(function.m_SmoothImageFwd), function.m_UnwarpImageFwd, function.m_CorrectedImageFwd );
  function.ComputeDeformedImage( x, -1, *(function.m_SmoothImageRev), function.m_UnwarpImageRev, function.m_CorrectedImageRev );

  // initialize gradient vector with derivative of image differences
  function.MakeGradientImage( x, +1, *(function.m_SmoothImageFwd), function.m_GradientImageFwd );
  function.MakeGradientImage( x, -1, *(function.m_SmoothImageRev), function.m_GradientImageRev );

  DataGrid::RegionType insideRegion = wholeImageRegion;  
  insideRegion.From()[function.m_PhaseEncodeDirection] += 1;
  insideRegion.To()[function.m_PhaseEncodeDirection] -= 1;
  size_t insideRegionSize = insideRegion.Size();

  const size_t nPixels = function.m_ImageGrid->GetNumberOfPixels();

  ap::real_value_type msd = 0;
  for ( RegionIndexIterator<DataGrid::RegionType> it( insideRegion ); it != it.end(); ++it )
    {
    // difference term derivative
    DataGrid::IndexType idx = it.Index();
    size_t px = sourceImage.GetOffsetFromIndex( idx );

    const Types::Coordinate diff = function.m_CorrectedImageFwd[px] - function.m_CorrectedImageRev[px];
    msd += MathUtil::Square( diff );
    g(1+px) = 2.0 * diff * (function.m_GradientImageFwd[px] + function.m_GradientImageRev[px]) / insideRegionSize;

    // add gradient terms for Jacobians
    //2*(j1(u, umm)*I1(x+um)-j2(u, umm)*I2(x-um)))*((diff(j1(u, umm), u))*I1(x+um)-(diff(j2(u, umm), u))*I2(x-um)
    idx[function.m_PhaseEncodeDirection] -= 1;
    px = sourceImage.GetOffsetFromIndex( idx );
    g(1+px) += 0.5 / insideRegionSize * (function.m_CorrectedImageFwd[px] - function.m_CorrectedImageRev[px]) * ( function.m_UnwarpImageFwd[px] + function.m_UnwarpImageRev[px] ); // "+" because partial J deriv is negative for Rev										       
    //(2*(j1(upp, u)*I1(x+up)-j2(upp, u)*I2(x-up)))*((diff(j1(upp, u), u))*I1(x+up)-(diff(j2(upp, u), u))*I2(x-up))
    idx[function.m_PhaseEncodeDirection] += 2;
    px = sourceImage.GetOffsetFromIndex( idx );
    // subtract second part because derivatives of J1 and J2 are negative here
    g(1+px) -= 0.5 / insideRegionSize * (function.m_CorrectedImageFwd[px] - function.m_CorrectedImageRev[px]) * ( function.m_UnwarpImageFwd[px] + function.m_UnwarpImageRev[px] ); // "+" because partial J deriv is negative for Rev
    }
  f = (msd /= insideRegionSize);

  // smoothness constraint and its derivative
  // compute smoothness term
  ap::real_value_type smooth = 0;
  const ap::real_value_type lambda2 = function.m_SmoothnessConstraintWeight;
  if ( lambda2 > 0 )
    {    
    insideRegion = wholeImageRegion;  
    insideRegion.From().AddScalar( 1 );
    insideRegion.To().AddScalar( -1 );
    insideRegionSize = insideRegion.Size();
    
    for ( RegionIndexIterator<DataGrid::RegionType> it( insideRegion ); it != it.end(); ++it )
      {
      const size_t ofs = 1 + sourceImage.GetOffsetFromIndex( it.Index() );
      for ( int dim = 0; dim < 3; ++dim )
	{
	const ap::real_value_type diff = x( ofs + sourceImage.m_GridIncrements[dim] ) - x( ofs - sourceImage.m_GridIncrements[dim] ) / 2;
	// increment smoothness term
	smooth += MathUtil::Square( diff );
	// increment relevant gradient elements
	g( ofs + sourceImage.m_GridIncrements[dim] ) += lambda2 * diff / insideRegionSize;
	g( ofs - sourceImage.m_GridIncrements[dim] ) -= lambda2 * diff / insideRegionSize;
	}
      }
    
    f += lambda2 * (smooth /= insideRegionSize);
    }
  
  DebugOutput( 5 ) << "f " << f << " msd " << msd << " smooth " << smooth << "\n";

}
