/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkScalarImageSimilarity.h"

#include <Base/cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{


ScalarImageSimilarity::ReturnType
ScalarImageSimilarity::GetRegionalMutualInformation
( const ScalarImage* image0, const ScalarImage* image1, const int radius )
{
  // printf("REGIONAL MUTUAL INFORMATION\n");
  const ScalarImage::IndexType& dims = image0->GetDims();
  int nrows = dims[1];
  int ncols = dims[0];
  
  int d = 2*radius+1;
  int N = (ncols-2*radius)*(nrows-2*radius);
  int d2 = d*d;
  int dim = 2*d2;
  
  Types::DataItem* pts = Memory::AllocateArray<Types::DataItem>( N*dim );
  
  
  Types::DataItem tmp1, tmp2;
  
  const TypedArray *data1 = image0->GetPixelData();
  const TypedArray *data2 = image1->GetPixelData();
  
  
  int idx,nidx,lidx;
  nidx = 0; idx = 0; lidx = 0;
  
// printf("creating neighborhood...\n");
// CREATE NEIGHBORHOOD VECTORS
  for (int i=radius; i<nrows-radius; i++) 
    {
    for (int j=radius; j<ncols-radius; j++) 
      {
      lidx = 0;
      for (int ii=-radius; ii<=radius; ii++) 
	{
	for (int jj=-radius; jj<=radius; jj++) 
	  {
	  idx = (i+ii)*ncols + j+jj;
	  
	  data1->Get(tmp1,idx);
	  data2->Get(tmp2,idx);
	  pts[lidx*N + nidx] = tmp1;
	  pts[(lidx+d2)*N + nidx] = tmp2;
	  lidx++;
	  }
        }
      nidx++;
      }
    }
// printf("done.\n");
  
  
// SUBTRACT MEAN
// printf("subtracting mean...\n");
  std::vector<double> mean(dim,0.0);
  
  for (int i=0; i<dim; i++) 
    {
    for (int j=0; j<N; j++) 
      {
      mean[i] += pts[i*N+j];
      }
    }
  for (int i=0; i<dim; i++) 
    {
    mean[i] /= N;
    }
  
  for (int i=0; i<dim; i++) 
    {
    for (int j=0; j<N; j++) 
      {
      pts[i*N+j] -= mean[i];
      }
    }
// printf("done.\n");
  
	

// CALCULATE JOINT COVARIANCE
// printf("calculating joint covariance...\n");
  Matrix2D<Types::DataItem> cM( dim, dim );
 
  ScalarImageSimilarity::ReturnType sum;
  int iN, jN;
  
  for (int i=0; i<dim; i++) 
    {
    for (int j=i; j<dim; j++) 
      {
      sum = 0.0;
      iN = i*N;
      jN = j*N;
      
      for (int k=0; k<N; k++) 
	{
	sum += pts[iN+k]*pts[jN+k];
	}
      cM[i][j] = sum/N;
      cM[j][i] = cM[i][j];
      // printf("cM[%d][%d] = %f\n", i,j,cM[i][j]);
      }
    }
  double dt3 = MathUtil::CholeskyDeterminant(cM, dim);
// printf("done.\n");
  
  
  
// CALCULATE MARGINAL COVARIANCES
// printf("calculating marginal covariances...\n");
  Matrix2D<Types::DataItem> mcM( dim/2, dim/2 );
 
  // image0's bloc
  for (int i=0; i<dim/2; i++) 
    {
    for (int j=0; j<dim/2; j++) 
      {
      mcM[i][j] = cM[i][j];
      }
    }
  double dt1 = MathUtil::CholeskyDeterminant(mcM, dim/2);
  
  // image1's bloc
  for (int i=0; i<dim/2; i++) 
    {
    for (int j=0; j<dim/2; j++) 
      {
      mcM[i][j] = cM[i+dim/2][j+dim/2];
      }
    }
  double dt2 =  MathUtil::CholeskyDeterminant(mcM, dim/2);
  
  double alpha = 1.41893853320467;
    
  double v1 = dim/2*alpha + .5*log(dt1);
  double v2 = dim/2*alpha + .5*log(dt2);
  double v3 = dim*alpha   + .5*log(dt3);
    
  // printf("v1=%f, v2=%f, v3=%f\n",v1,v2,v3);
    
  return ((ScalarImageSimilarity::ReturnType) (v1+v2-v3));
}

} // namespace cmtk
