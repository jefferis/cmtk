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

#include <cmtkTypedArraySimilarity.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{


TypedArraySimilarity::ReturnType
TypedArraySimilarity::GetMutualInformation
( const std::vector<const TypedArray*>& data0, 
  const std::vector<const TypedArray*>& data1,
  const bool normalized )
{
  const size_t N = data0[0]->GetDataSize();
  const size_t dim0 = data0.size();
  const size_t dim1 = data1.size();
  const size_t dim = dim0 + dim1;
  
  double* pts = Memory::AllocateArray<double>( N*dim );  
  
  Types::DataItem tmp;
  
// CREATE NEIGHBORHOOD VECTORS
  for ( size_t nidx = 0; nidx < N; ++nidx )
    {
    for ( size_t lidx = 0; lidx < dim0; ++lidx )
      {
      data0[lidx]->Get(tmp,nidx);
      pts[lidx * N + nidx] = tmp;
      }
    for ( size_t lidx = 0; lidx < dim1; ++lidx )
      {
      data1[lidx]->Get(tmp,nidx);
      pts[(dim0 + lidx) * N + nidx] = tmp;
      }    
    }
  
  
// SUBTRACT MEAN
  std::vector<double> mean(dim,0.0);
  
  for (size_t i=0; i<dim; i++) 
    {
    for (size_t j=0; j<N; j++) 
      {
      mean[i] += pts[i*N+j];
      }
    }
  for (size_t i=0; i<dim; i++) 
    {
    mean[i] /= N;
    }
  
  for (size_t i=0; i<dim; i++) 
    {
    for (size_t j=0; j<N; j++) 
      {
      pts[i*N+j] -= mean[i];
      }
    }

// CALCULATE JOINT COVARIANCE
// printf("calculating joint covariance...\n");
  Matrix2D<Types::DataItem> cM( dim, dim );
 
  double sum;
  size_t iN, jN;
  
  for (size_t i=0; i<dim; i++) 
    {
    for (size_t j=i; j<dim; j++) 
      {
      sum = 0.0;
      iN = i*N;
      jN = j*N;
      
      for (size_t k=0; k<N; k++) 
	{
	sum += pts[iN+k]*pts[jN+k];
	}
      cM[i][j] = sum/N;
      cM[j][i] = cM[i][j];
      }
    }
  const double dt3 =  MathUtil::CholeskyDeterminant(cM, dim);
  
// CALCULATE MARGINAL COVARIANCES
// printf("calculating marginal covariances...\n");
  Matrix2D<Types::DataItem> mcM0( dim0, dim0 );
  
  // image0's bloc
  for (size_t i=0; i<dim0; i++) 
    {
    for (size_t j=0; j<dim0; j++) 
      {
      mcM0[i][j] = cM[i][j];
      }
    }

  const double dt1 = MathUtil::CholeskyDeterminant(mcM0, dim0);

    // image1's bloc
  Matrix2D<Types::DataItem> mcM1( dim1, dim1 );
  
  // image0's bloc
  for (size_t i=0; i<dim1; i++) 
    {
    for (size_t j=0; j<dim1; j++) 
      {
      mcM1[i][j] = cM[i+dim0][j+dim0];
      }
    }

  const double dt2 = MathUtil::CholeskyDeterminant(mcM1, dim1);

  const double alpha = 1.41893853320467;
    
  const double v1 = dim0*alpha + .5*log(dt1);
  const double v2 = dim1*alpha + .5*log(dt2);
  const double v3 = dim *alpha + .5*log(dt3);
    
  if ( normalized )
    return (v1+v2)/v3;
  else
    return v1+v2-v3;
}

} // namespace cmtk
