/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <algorithm>

#include <Base/cmtkHistogramBase.h>

#include <System/cmtkConsole.h>

#ifdef CMTK_USE_GCD
#  include <dispatch/dispatch.h>
#endif

template<class TPSF>
void
DeblurringVolumeReconstruction<TPSF>
::Blur( const ap::real_1d_array& reconstructedPixelArray )
{
  this->m_InterpolatedPassImages.clear();
  
  const UniformVolume* correctedImage = this->m_CorrectedImage;
  const DataGrid::IndexType& correctedImageDims = correctedImage->GetDims();

  for ( int pass = 0; pass < this->m_NumberOfPasses; ++pass )
    {
    const UniformVolume* passImage = this->m_OriginalPassImages[pass];
    UniformVolume::SmartPtr result( passImage->CloneGrid() );
    result->CreateDataArray( TYPE_FLOAT, true/*setToZero*/ );
    
    const DataGrid::IndexType& passImageDims = passImage->GetDims();
    const int passImageDimsX = passImageDims[0], passImageDimsY = passImageDims[1], passImageDimsZ = passImageDims[2];
    const size_t passImageDimsXY = passImageDimsX*passImageDimsY;
    const size_t passImageDimsXYZ = passImageDimsXY*passImageDimsZ;

    const TPSF* passImagePSF = this->m_PassImagePSF[pass];
    const AffineXform* correctedToPassXform = AffineXform::SmartPtr::DynamicCastFrom( this->m_TransformationsToPassImages[pass] );
    const AffineXform* passToCorrectedXform = AffineXform::SmartPtr::DynamicCastFrom( this->m_TransformationsToPassImages[pass] )->GetInverse();

#ifdef CMTK_USE_GCD
    const size_t stride = passImageDimsXYZ / (2 * Threads::GetNumberOfProcessors() );
    dispatch_apply( passImageDimsXYZ / stride, , dispatch_get_global_queue(0, 0), ^(size_t i ) 
		    { const size_t last = std::min( stride * (i+1), passImageDimsXYZ ); for ( size_t offset = i * stride; offset < last; ++offset )
#else
#pragma omp parallel for
    for ( size_t offset = 0; offset < passImageDimsXYZ; ++offset )
#endif
      {
      const int curPassVox[3] = { (offset % passImageDimsXY) % passImageDimsX, (offset % passImageDimsXY) / passImageDimsX, (offset / passImageDimsXY) };
      
      Types::DataItem blurredData = 0;
      Types::DataItem totalWeight = 0;
      
      /* Get a bounding box for the transformed neighborhood (pass- to corrected-image)
           */
      const UniformVolume::CoordinateVectorType curPassVec3D = passImage->GetGridLocation( curPassVox[0], curPassVox[1], curPassVox[2] );
      int corrBoundingBox[6];
      this->GetBoundingBoxOfXformedPassNeighborhood( corrBoundingBox, correctedImage, curPassVec3D, passImagePSF, passToCorrectedXform, correctedImageDims );
      
      /* Iterate through the bounding box of the blur-weights matrix,
           * assembling the weighted sum of the neighbors, and the sum of weights
           */
      Types::DataItem data;
      for (int k = corrBoundingBox[2]; k <= corrBoundingBox[5]; ++k)
	{
	for (int j = corrBoundingBox[1]; j <= corrBoundingBox[4]; ++j)
	  {
	  for (int i = corrBoundingBox[0]; i <= corrBoundingBox[3]; ++i)
	    {
	    UniformVolume::CoordinateVectorType curNeighbVec3D;
	    Types::Coordinate from[3], to[3];
	    curNeighbVec3D = correctedImage->GetGridLocation( i, j, k );
	    correctedToPassXform->ApplyInPlace( curNeighbVec3D );
            
	    int neighborPassVox[3];
	    if ( passImage->FindVoxel( curNeighbVec3D, neighborPassVox, from, to ) )
	      {
	      Types::Coordinate weightIJK = 1;
	      for ( int dim = 0; dim < 3; ++dim )
		{
		const Types::Coordinate displacement = curNeighbVec3D[dim] - curPassVec3D[dim];
		weightIJK *= passImagePSF->GetWeight( dim, displacement );
		}
	      
	      totalWeight += weightIJK;
	      data = reconstructedPixelArray( 1+correctedImage->GetOffsetFromIndex( i, j, k ) );                 
	      blurredData += data * weightIJK;
	      }
	    }
	  }
	}
      /* Assign the blurred value to the result */
      if ( totalWeight == 0 )
	result->GetData()->SetPaddingAt( offset );
      else
	result->SetDataAt( blurredData / totalWeight, offset );
      }
    this->m_InterpolatedPassImages.push_back( result );
    }
#ifdef CMTK_USE_GCD
		    } );
#endif
}

template <class TPSF>
void
DeblurringVolumeReconstruction<TPSF>
::ComputeErrorGradientImage( ap::real_1d_array& g )
{
  const UniformVolume* correctedImage = this->m_CorrectedImage;
  const size_t numberOfPixels = correctedImage->GetNumberOfPixels();
  for ( size_t i = 1; i <= numberOfPixels; ++i )
    g(i) = 0;

  const DataGrid::IndexType& correctedImageDims = correctedImage->GetDims();
  const int correctedImageDimsX = correctedImageDims[0], correctedImageDimsY = correctedImageDims[1];
  const int correctedImageDimsXY = correctedImageDimsX*correctedImageDimsY;

  for ( int pass = 0; pass < this->m_NumberOfPasses; ++pass )
    {
    const Types::Coordinate passImageWeight = this->m_PassWeights[pass];      
    if ( passImageWeight > 0 )
      {
      const UniformVolume* differencePassImage = this->m_DifferencePassImages[pass];
      const UniformVolume* blurredPassImage = this->m_InterpolatedPassImages[pass];
      
      const TPSF* passImagePSF = this->m_PassImagePSF[pass];
      const AffineXform* transformationToPassImage = AffineXform::SmartPtr::DynamicCastFrom( this->m_TransformationsToPassImages[pass] );
      
#pragma omp parallel for
      for ( size_t offset = 0; offset < numberOfPixels; ++offset )
	{
	const int corrCenterVox[3] = { (offset % correctedImageDimsXY) % correctedImageDimsX, (offset % correctedImageDimsXY) / correctedImageDimsX, (offset / correctedImageDimsXY) };

        UniformVolume::CoordinateVectorType curCenterVec3D = correctedImage->GetGridLocation( corrCenterVox[0], corrCenterVox[1], corrCenterVox[2] );
        transformationToPassImage->ApplyInPlace( curCenterVec3D );

	/* compute neighborhood in blurred (pass) image from which blurred pixels 
           * are affected by current corrected image pixel
           */
	int passRegionCorners[2][3];
	for ( int dim = 0; dim < 3; ++dim )
	  {
	  passRegionCorners[0][dim] = blurredPassImage->GetCoordIndex( dim, curCenterVec3D[dim] - passImagePSF->GetTruncationRadius( dim ) );
	  passRegionCorners[1][dim] = blurredPassImage->GetCoordIndex( dim, curCenterVec3D[dim] + passImagePSF->GetTruncationRadius( dim ) );
	  }
	
	/* Iterate through the neighborhood in the blurred (pass) image,
           * assembling a value for the gradient at this offset of the corrected
           * image
           */
	Types::DataItem gradientData = 0;
	for (int k = passRegionCorners[0][2]; k <= passRegionCorners[1][2]; ++k)
	  {
	  for (int j = passRegionCorners[0][1]; j <= passRegionCorners[1][1]; ++j)
	    {
	    for (int i = passRegionCorners[0][0]; i <= passRegionCorners[1][0]; ++i)
	      {
	      Types::DataItem differenceData;
	      // if differenceData==0, we can skip this, too
	      if ( differencePassImage->GetDataAt( differenceData, i, j, k ) && (differenceData !=0) ) 
		{
		Types::DataItem weight = 0;
		
		const UniformVolume::CoordinateVectorType passNeighbVec3D = blurredPassImage->GetGridLocation( i, j, k );
                
		weight = passImageWeight;
		for ( int dim = 0; dim < 3; ++dim )
		  {
		  const Types::Coordinate displacement = curCenterVec3D[dim] - passNeighbVec3D[dim];
		  weight *= passImagePSF->GetWeight( dim, displacement );
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

template<class TPSF>
void
DeblurringVolumeReconstruction<TPSF>
::GetBoundingBoxOfXformedPassNeighborhood
( int* correctedImageBoundingBox, const UniformVolume* correctedImage, const Vector3D& currentPassVoxel, const TPSF* psf,
  const AffineXform* passToCorrectedXform, const DataGrid::IndexType& correctedImageDims ) const
{
  /* Compute the blurring neighborhood in the pass image 
   */
  Types::Coordinate corners[2][3];
  for ( int dim = 0; dim < 3; ++dim )
    {
    const Types::Coordinate radius = psf->GetTruncationRadius( dim );
    corners[0][dim] = currentPassVoxel[dim] - radius;
    corners[1][dim] = currentPassVoxel[dim] + radius;
    }

  /* Make corner vectors out of that 2-D array and put them in a 1-D array,
   * transforming each to the space of the corrected image.
   */
  Vector3D corners3D[8];
  int neighborIdx = 0;
  for ( int a = 0; a < 2 ; ++a )
    {
    for ( int b = 0; b < 2 ; ++b )
      {
      for ( int c = 0; c < 2; ++c, ++neighborIdx )
        {
        corners3D[neighborIdx][0] = corners[a][0];
        corners3D[neighborIdx][1] = corners[b][1];
        corners3D[neighborIdx][2] = corners[c][2];
	passToCorrectedXform->ApplyInPlace( corners3D[neighborIdx] );
        }
      }
    }

  /* Compute the bounding box of that transformed region 
   */
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

  /* Put the voxel indices corresponding to the corners of that bounding box 
   * into correctedImageBoundingBox
   */
  correctedImage->GetVoxelIndexNoBounds( bboxMin, correctedImageBoundingBox+0 );
  correctedImage->GetVoxelIndexNoBounds( bboxMax, correctedImageBoundingBox+3 );

  /* Clip bounding box (now correctedImageBoundingBox) against corrected image boundaries 
   */
  for ( int dim = 0; dim < 3; ++dim )
    {
    correctedImageBoundingBox[dim] = std::max( correctedImageBoundingBox[dim], 0 );
    // increment upper indexes by one to compensate for floating point truncation in pixel index lookup.
    correctedImageBoundingBox[dim+3] = std::min( correctedImageBoundingBox[dim+3]+1, correctedImageDims[dim]-1 );
    }
}

template<class TPSF>
void
DeblurringVolumeReconstruction<TPSF>
::FunctionAndGradient
::Evaluate( const ap::real_1d_array& x, ap::real_value_type& f, ap::real_1d_array& g )
{
  this->m_Function->Blur( x );
  this->m_Function->ComputeApproximationError();
  this->m_Function->ComputeErrorGradientImage( g );
  
  const ap::real_value_type msd = f = this->m_Function->GetMeanSquaredError();

  ap::real_value_type lnorm = 0;
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
  
  StdOout << "f " << f << " MSD " << msd
	  << " MAX " << this->m_Function->GetMaximumError() 
	  << " KLD " << this->m_Function->GetOriginalToCorrectedImageKLD( x )
	  << " LNORM " << lnorm << "\n";
}
