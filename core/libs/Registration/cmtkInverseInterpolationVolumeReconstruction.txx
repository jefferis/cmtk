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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <algorithm>

#include <cmtkHistogramBase.h>

namespace
cmtk 
{

/** \addtogroup Registration */
//@{

template <class TInterpolator>
void
InverseInterpolationVolumeReconstruction<TInterpolator>
::Interpolation( const ap::real_1d_array& reconstructedPixelArray )
{
  this->m_InterpolatedPassImages.clear();
  for ( int pass = 0; pass < this->m_NumberOfPasses; ++pass )
    {
    const UniformVolume* passImage = this->m_OriginalPassImages[pass];
    UniformVolume::SmartPtr result( passImage->CloneGrid() );
    result->CreateDataArray( TYPE_FLOAT, true/*setToZero*/ );
    
    const int *passImageDims = passImage->GetDims();
    const int passImageDimsX = passImageDims[0], passImageDimsY = passImageDims[1], passImageDimsZ = passImageDims[2];
    
    const int *correctedImageDims = this->m_CorrectedImage->GetDims();
    const int correctedImageDimsX = correctedImageDims[0], correctedImageDimsY = correctedImageDims[1], correctedImageDimsZ = correctedImageDims[2];    

    const AffineXform* passXform = dynamic_cast<const AffineXform*>( this->m_TransformationsToPassImages[pass].GetPtr() )->GetInverse();
    
#pragma omp parallel for
    for ( int x = 0; x < passImageDimsX; ++x )
      {
      Vector3D v;
      int correctedImageGridPoint[3];
      Types::Coordinate from[3], to[3];
    
      for ( int y = 0; y < passImageDimsY; ++y )
	{
        for ( int z = 0; z < passImageDimsZ; ++z )
	  {
          Types::DataItem interpolatedData = 0;
	  Types::Coordinate totalWeight = 0;
          passImage->GetGridLocation( v, x, y, z );
	  passXform->ApplyInPlace( v );

          if ( this->m_CorrectedImage->FindVoxel( v, correctedImageGridPoint, from, to ) )
	    {
            const int xx = correctedImageGridPoint[0] + 1 - TInterpolator::RegionSizeLeftRight;
            const int yy = correctedImageGridPoint[1] + 1 - TInterpolator::RegionSizeLeftRight;
            const int zz = correctedImageGridPoint[2] + 1 - TInterpolator::RegionSizeLeftRight;

            Types::DataItem data;
            Types::Coordinate difference[3];
  	    Types::Coordinate interpolationWeights[3][2 * TInterpolator::RegionSizeLeftRight];
	    for ( int n = 0; n < 3; ++n )
	      {
	      difference[n] = (v[n] - from[n]) / (to[n] - from[n]);
	      for ( int m = 1-TInterpolator::RegionSizeLeftRight; m <= TInterpolator::RegionSizeLeftRight; ++m )
		{
		interpolationWeights[n][m+TInterpolator::RegionSizeLeftRight-1] = TInterpolator::GetWeight(m, difference[n]);
		}
	      }
	    
	    const int iMin = std::max( 0, -xx );
	    const int iMax = std::min( 2 * TInterpolator::RegionSizeLeftRight, correctedImageDimsX - xx );

	    const int jMin = std::max( 0, -yy );
	    const int jMax = std::min( 2 * TInterpolator::RegionSizeLeftRight, correctedImageDimsY - yy );

	    const int kMin = std::max( 0, -zz );
	    const int kMax = std::min( 2 * TInterpolator::RegionSizeLeftRight, correctedImageDimsZ - zz );

	    for ( int k = kMin; k < kMax; ++k )
	      {
	      for ( int j = jMin; j < jMax; ++j )
		{
		const Types::Coordinate weightJK = interpolationWeights[1][j] * interpolationWeights[2][k];
		for ( int i = iMin; i < iMax; ++i )
		  {
		  const Types::Coordinate weightIJK = interpolationWeights[0][i] * weightJK;
		  data = reconstructedPixelArray( 1+this->m_CorrectedImage->GetOffsetFromIndex( xx + i, yy + j, zz + k ) );

		  interpolatedData = interpolatedData + data * weightIJK;
		  totalWeight += weightIJK;
		  }
		}
	      }
	    }
	  
	  if ( totalWeight == 0 )
	    result->GetData()->SetPaddingAt( result->GetOffsetFromIndex( x, y, z ) );
	  else
	    result->SetDataAt( interpolatedData / totalWeight, x, y, z );
	  }
	}
      }
    this->m_InterpolatedPassImages.push_back( result );
    }
}

template <class TInterpolator>
void
InverseInterpolationVolumeReconstruction<TInterpolator>
::ComputeErrorGradientImage( ap::real_1d_array& g )
{
  const UniformVolume* correctedImage = this->m_CorrectedImage;
  const size_t numberOfPixels = correctedImage->GetNumberOfPixels();
  for ( size_t i = 1; i <= numberOfPixels; ++i )
    g(i) = 0;

  const int *correctedImageDims = correctedImage->GetDims();
  const int correctedImageDimsX = correctedImageDims[0], correctedImageDimsY = correctedImageDims[1];
  const int correctedImageDimsXY = correctedImageDimsX*correctedImageDimsY;

  for ( int pass = 0; pass < this->m_NumberOfPasses; ++pass )
    {
    const UniformVolume* differencePassImage = this->m_DifferencePassImages[pass];
    const UniformVolume* interpolatedPassImage = this->m_InterpolatedPassImages[pass];

    const AffineXform* transformationToPassImage = dynamic_cast<const AffineXform*>( this->m_TransformationsToPassImages[pass].GetPtr() );
    const AffineXform* inverseTransformationToPassImage = transformationToPassImage->GetInverse();
    const int *passImageDims = this->m_InterpolatedPassImages[pass]->GetDims();

    const Types::Coordinate passImageWeight = this->m_PassWeights[pass];

    if ( passImageWeight > 0 )
      {
#pragma omp parallel for
      for ( size_t offset = 0; offset < numberOfPixels; ++offset )
	{
	const int correctedImageCurrentGridPoint[3] = 
	  { (offset % correctedImageDimsXY) % correctedImageDimsX, (offset % correctedImageDimsXY) / correctedImageDimsX, (offset / correctedImageDimsXY) };
	
	int passImageDependentRegion[6];
	this->GetPassImageDependentPixelRegion
	  ( passImageDependentRegion, correctedImage, correctedImageCurrentGridPoint, interpolatedPassImage, transformationToPassImage, passImageDims );

	Types::DataItem gradientData = 0;
	Types::Coordinate from[3], to[3];
	for (int k = passImageDependentRegion[2]; k <= passImageDependentRegion[5]; ++k)
	  {
	  for (int j = passImageDependentRegion[1]; j <= passImageDependentRegion[4]; ++j)
	    {
	    for (int i = passImageDependentRegion[0]; i <= passImageDependentRegion[3]; ++i)
	      {
	      Types::DataItem differenceData;
	      if ( differencePassImage->GetDataAt( differenceData, i, j, k ) && (differenceData !=0) ) // if differenceData==0, we can skip this, too
		{
		Types::DataItem weight = 0;
	      
		Vector3D v;
		interpolatedPassImage->GetGridLocation( v, i, j, k );
		inverseTransformationToPassImage->ApplyInPlace( v );
	      
		int correctedImageGridPoint[3];
		if ( correctedImage->FindVoxel( v, correctedImageGridPoint, from, to ) )
		  {
		  const int correctedImageDifference[3] = 
		    {
		      correctedImageCurrentGridPoint[0] - correctedImageGridPoint[0],
		      correctedImageCurrentGridPoint[1] - correctedImageGridPoint[1],
		      correctedImageCurrentGridPoint[2] - correctedImageGridPoint[2],
		    };
		
		  if ( correctedImageDifference[0] > -TInterpolator::RegionSizeLeftRight && correctedImageDifference[1] > -TInterpolator::RegionSizeLeftRight && 
		       correctedImageDifference[2] > -TInterpolator::RegionSizeLeftRight && correctedImageDifference[0] <= TInterpolator::RegionSizeLeftRight &&
		       correctedImageDifference[1] <= TInterpolator::RegionSizeLeftRight && correctedImageDifference[2] <= TInterpolator::RegionSizeLeftRight )
		    {
		    weight = passImageWeight;
		    for ( int n = 0; n < 3; ++n )
		      {
		      const Types::Coordinate relative = (v[n] - from[n]) / (to[n] - from[n]);
		      weight *= TInterpolator::GetWeight( correctedImageDifference[n], relative );
		      }
		    }
		  }
	      
		gradientData += weight * differenceData;
		}
	      }
	    }
	  }
	
	if ( this->m_FourthOrderError )
	  g(offset+1) += 4 * (gradientData*gradientData*gradientData) / numberOfPixels;
	else
	  g(offset+1) += 2 * gradientData / numberOfPixels;
	}
      }
    }
}

template<class TInterpolator>
void
InverseInterpolationVolumeReconstruction<TInterpolator>
::GetPassImageDependentPixelRegion
( int* region, const UniformVolume* correctedImage, const int* currentCorrectedGridPoint, 
  const UniformVolume* passImage, const AffineXform* transformationToPassImage, const int* passImageDims )
{
  // compute neighborhood in corrected image from which interpolated pixels are affected by current corrected image pixel
  int corners[2][3];
  for ( int dim = 0; dim < 3; ++dim )
    {
    corners[0][dim] = std::max( currentCorrectedGridPoint[dim] - TInterpolator::RegionSizeLeftRight, 0 );
    corners[1][dim] = std::min( currentCorrectedGridPoint[dim] + TInterpolator::RegionSizeLeftRight, correctedImage->GetDims(dim)-1 );
    }

  Vector3D corners3D[8];
  int neighborIdx = 0;
  for ( int a = 0; a < 2 ; ++a )
    {
    for ( int b = 0; b < 2 ; ++b )
      {
      for ( int c = 0; c < 2; ++c, ++neighborIdx )
	{
	correctedImage->GetGridLocation( corners3D[neighborIdx], corners[a][0], corners[b][1], corners[c][2] );
	transformationToPassImage->ApplyInPlace( corners3D[neighborIdx] );
	}
      }
    }
  
  // now get bounding box of transformed region that is aligned with pass image
  Vector3D bboxMin = corners3D[0];
  Vector3D bboxMax = corners3D[0];
  for ( int m = 1; m < 8; ++m )
    {
    for ( int o = 0; o < 3; ++o )
      {
      if ( corners3D[m][o] < bboxMin[o] )
	{
	bboxMin[o] = corners3D[m][o];
	}
      else
	if ( corners3D[m][o] > bboxMax[o] )
	  {
	  bboxMax[o] = corners3D[m][o];
	  }
      }
    }

  // get pass image pixels indexes corresponding to bounding box
  passImage->GetVoxelIndexNoBounds( bboxMin, region+0 );
  passImage->GetVoxelIndexNoBounds( bboxMax, region+3 );
  
  // clip bounding box against pass image boundaries
  for ( int dim = 0; dim < 3; ++dim )
    {
    region[dim] = std::max( region[dim], 0 );
    // increment upper indexes by one to compensate for floating point truncation in pixel index lookup.
    region[dim+3] = std::min( region[dim+3]+1, passImageDims[dim]-1 );
    }
}

template<class TInterpolator>
void
InverseInterpolationVolumeReconstruction<TInterpolator>
::FunctionAndGradient
::Evaluate( const ap::real_1d_array& x, double& f, ap::real_1d_array& g )
{
  this->m_Function->Interpolation( x );
  this->m_Function->ComputeApproximationError();
  this->m_Function->ComputeErrorGradientImage( g );
  const double msd = f = this->m_Function->GetMeanSquaredError();

  double lnorm = 0;
  if ( this->m_Function->m_ConstraintWeightLNorm > 0 )
    {
    f += this->m_Function->m_ConstraintWeightLNorm * (lnorm = this->m_Function->ComputeCorrectedImageLaplacianNorm( x ));
    this->m_Function->AddLaplacianGradientImage( g, x, this->m_Function->m_ConstraintWeightLNorm );
    }

  if ( this->m_Function->GetMaximumError() <= this->m_Function->m_LowestMaxError )
    {
    this->m_Function->m_LowestMaxError = this->m_Function->GetMaximumError();
    const int numberOfPixels = this->m_Function->m_CorrectedImage->GetNumberOfPixels();
    for ( int i = 1; i <= numberOfPixels; ++i )
      this->m_Function->m_CorrectedImage->SetDataAt( x(i), i-1 );
    this->m_Function->m_LowestMaxErrorImage = UniformVolume::SmartPtr( this->m_Function->m_CorrectedImage->Clone( true /*copyData*/ ) );
    }
  
  std::cerr << "f " << f << " MSD " << msd
	    << " MAX " << this->m_Function->GetMaximumError() 
	    << " KLD " << this->m_Function->GetOriginalToCorrectedImageKLD( x )
	    << " LNORM " << lnorm << std::endl;
}

} // namespace cmtk
