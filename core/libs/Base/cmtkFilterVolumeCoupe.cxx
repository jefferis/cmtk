/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <cmtkFilterVolume.h>

namespace
cmtk
{

Types::DataItem
FilterVolume::Mean
( CoupeBlock items )
{
  Types::DataItem sum = 0.0;
  for ( int i = 0; i < COUPE_BLOCK_SIZE; i++ )
    sum += items[i];
  return sum / COUPE_BLOCK_SIZE;
}

Types::DataItem
FilterVolume::Variance
( CoupeBlock items, const Types::DataItem mean )
{
  Types::DataItem sum = 0.0;
  for ( int i = 0; i < COUPE_BLOCK_SIZE; i++ )
    sum += pow( items[i] - mean, 2 );
  return sum / COUPE_BLOCK_SIZE;
}

void
FilterVolume::BlockAddInPlace
( CoupeBlock v1, CoupeBlock v2 )
{
  
  for ( int i = 0; i < COUPE_BLOCK_SIZE; i++ )
    v1[i] += v2[i] ;
}

void
FilterVolume::BlockSubtract
( CoupeBlock diff, CoupeBlock v1, CoupeBlock v2 )
{
  for ( int i = 0; i < COUPE_BLOCK_SIZE; i++ )
    diff[i] = v1[i] - v2[i];
}

void
FilterVolume::BlockConstMult
( CoupeBlock prod, CoupeBlock items, const Types::DataItem mult )
{
  for ( int i = 0; i < COUPE_BLOCK_SIZE; i++ )
    prod[i] = items[i] * mult;
}

double
FilterVolume::BlockSquaredDistance
( CoupeBlock centerBlock, CoupeBlock outerBlock )
{
  CoupeBlock diff; 
  BlockSubtract( diff, centerBlock, outerBlock );   
  
  Types::DataItem sum = 0.0;
  for ( int i = 0; i < COUPE_BLOCK_SIZE; i++ )
    sum += pow( diff[i], 2 );
  return ( sum );  // Coupe paper uses squared distance
}

void
FilterVolume::GetCoupeBlock
( CoupeBlock block,
  const TypedArray* data, const int* dims,
  const int x, const int y, const int z )
{

  const int blockRadius = COUPE_BLOCK_RADIUS;
  int curBlockSlot = 0;
  /*  No bounds-checking here.  The program should fail if 
   *  a block is requested too close to the edge of the image.
   */
  for ( int k = z - blockRadius; k <= z + blockRadius; ++k )
    for ( int j = y - blockRadius; j <= y + blockRadius; ++j )
      for ( int i = x - blockRadius; i <= x + blockRadius; ++i ) 
        {
        int offset = i + dims[AXIS_X] * ( j + dims[AXIS_Y] * k );
        Types::DataItem value;
        data->Get( value, offset );
        block[curBlockSlot++] = value;
        }
}

double 
FilterVolume::ComputeCoupeWeight
( const Types::DataItem smoothingParam, 
  CoupeBlock centerBlock, 
  CoupeBlock outerBlock )
{
  const double numerator = BlockSquaredDistance( centerBlock, outerBlock );
  
  double weight = 1.0;
  if ( numerator != 0 )
    weight = exp( -(log(numerator) / smoothingParam) );
  //std::cout << numerator << ", " << smoothingParam << std::endl;
  //std::cout << weight << "\n";
  return weight;
}

void
FilterVolume::ComputeNLWithinWindow
( CoupeBlock NL,
  const TypedArray* blockLocations,
  const TypedArray* data, const int* dims, const Types::DataItem smoothingParam,
  const int x, const int y, const int z, 
  const int windowRadius, 
  const float, // beta 
  const TypedArray* localMeansMap, 
  const TypedArray* localVariancesMap, 
  CoupeBlock centerBlock )
{

  /*  These two constants were determined experimentally
   *  by Coupe, et al.  (Section V-C, Fig. 8.)
   */
  const float Mu1 = 0.95;
  const float Sigma1SQ = 0.5;
  
  /*  In this loop, build a set of weights from only
   *  that portion of the neighborhood-window that falls
   *  within the image boundaries.
   */
  int curNeighbor = 0;
  int centerBlockPos = -1;
  const int windowSize = (2 * windowRadius + 1) * (2 * windowRadius + 1) * (2 * windowRadius + 1);
  const int blockRadius = COUPE_BLOCK_RADIUS;
  int centerOffset = x + dims[AXIS_X] * ( y + dims[AXIS_Y] * z );
  Types::DataItem centerBlockMean, centerBlockVariance;
  localMeansMap->Get( centerBlockMean, centerOffset );
  localVariancesMap->Get( centerBlockVariance, centerOffset );
  CoupeBlock curBlock;
  int offset = 0;
  Types::DataItem blockAtCurVox = 0.0;
#ifdef CMTK_VAR_AUTO_ARRAYSIZE
  Types::DataItem neighborBlocks[windowSize][COUPE_BLOCK_SIZE];
  Types::DataItem neighborWeights[windowSize];
#else
  Matrix2D<Types::DataItem> neighborBlocks( windowSize, COUPE_BLOCK_SIZE );
  std::vector<Types::DataItem> neighborWeights( windowSize );
#endif
  Types::DataItem curBlockMean, curBlockVariance;
  Types::DataItem ratioBlockMean, ratioBlockVariance;
  for ( int k = z - windowRadius; k <= z + windowRadius; k++ )
    if ( ( k >= blockRadius ) && ( k < dims[AXIS_Z] - blockRadius ) )
      {
      for ( int j = y - windowRadius; j <= y + windowRadius; j++ )
        if ( ( j >= blockRadius ) && ( j < dims[AXIS_Y] - blockRadius ) )
          {
          for ( int i = x - windowRadius; i <= x + windowRadius; i++ )
            if ( ( i >= blockRadius ) && ( i < dims[AXIS_X] - blockRadius ) )
              {
              offset = i + dims[AXIS_X] * ( j + dims[AXIS_Y] * k );
              blockLocations->Get( blockAtCurVox, offset ); 
              if ( blockAtCurVox ) // if there's a block here
                {
                GetCoupeBlock( curBlock, data, dims, i, j, k );
                localMeansMap->Get( curBlockMean, offset );
                localVariancesMap->Get( curBlockVariance, offset );
                ratioBlockMean = centerBlockMean / curBlockMean;
                /*  This block is to handle the 0:0 ratio without divide-by-zero
                 */
                if ( curBlockVariance != 0 ) 
                  {
                  ratioBlockVariance = centerBlockVariance / curBlockVariance;
                  }
                else
                  {
                  ratioBlockVariance = ( centerBlockVariance == 0 ) ? 1.0 : 0.0;
                  }

                if ( ( Mu1 <= ratioBlockMean ) && ( ratioBlockMean <= ( 1 / Mu1 ) )
                  && ( Sigma1SQ <= ratioBlockVariance ) && ( ratioBlockVariance <= ( 1 / Sigma1SQ ) ) )
                  {
                  /*  Handle blocks that aren't the center block here,
                   *  handle the center block in 'else' below
                   */
                  if ( ( ( i != x ) || ( j != y ) || ( k != z ) ) ) 
                    {
                    double weight = ComputeCoupeWeight( smoothingParam, centerBlock, curBlock );
                    if ( weight > 0 )
                      {
                      if ( ( i == 91 ) && ( j == 1 ) && ( k == 93 ) ) 
                        {
                        Types::DataItem tmp;
                        data->Get( tmp, offset );
                        std::cout << i << "," << j << "," << k << "\t";
                        std::cout << tmp << "\t" << curBlockMean << "\t" << weight << "\n";
                        }
                      memcpy( neighborBlocks[curNeighbor], curBlock, sizeof( curBlock ) );
                      neighborWeights[curNeighbor] = weight;
                      }
                    }
                  /*  Handle center block
                   */
                  else
                    {
                    if ( ( i == 91 ) && ( j == 1 ) && ( k == 93 ) ) 
                      {
                      Types::DataItem tmp;
                      data->Get( tmp, offset );
                      std::cout << i << "," << j << "," << k << "\t" << tmp << "\t";
                      std::cout << curBlockMean << "\t..." << "\n";
                      }
                    memcpy( neighborBlocks[curNeighbor], curBlock, sizeof( curBlock ) );
                    
                    /*  Store the position of the center block
                     *  w.r.t. the list of used neighbors, for use
                     *  in the weighting clause below.
                     */
                    centerBlockPos = curNeighbor;
                    }
                  curNeighbor++;
                  }
                }
              }
          }
      }
  int numNeighborsUsed = curNeighbor + 1;

  /*  As per Manjon et al, 2008, set the weight of the center
   *  block to equal the largest of the neighbor-weights.
   */
  if ( numNeighborsUsed < 2 )
    {
    neighborWeights[centerBlockPos] = 1.0;
    }
  else
    {
    Types::DataItem maxWt = 0.0;
    for ( int i = 0; i < numNeighborsUsed; i++ )
      maxWt = ( neighborWeights[i] > maxWt ) ? neighborWeights[i] : maxWt;
    //neighborWeights[centerBlockPos] = 1.0;
    neighborWeights[centerBlockPos] = maxWt;
    }

  /*  Sum the weights for normalization
   */
  Types::DataItem weightsSum = 0.0;
  for ( int i = 0; i < numNeighborsUsed; i++ )
    weightsSum += neighborWeights[i];

  for ( int i = 0; i < COUPE_BLOCK_SIZE; i++ )
    NL[i] = 0.0;

  //if ( weightsSum > 1e-9 ) 
  if ( weightsSum != 0 ) 
    {
    /*  Multiply each neighborBlock by its normalized weight,
     *  then add that to the output NL value.
     */  
    CoupeBlock weightedCurBlock;
    Types::DataItem tmp = 0.0;
    for ( int i = 0; i < numNeighborsUsed; i++ )
      {
      Types::DataItem normalizedWeight = 0;
        normalizedWeight = neighborWeights[i] / weightsSum;

      tmp += normalizedWeight;
      BlockConstMult( weightedCurBlock, neighborBlocks[i], normalizedWeight );
      BlockAddInPlace( NL, weightedCurBlock );
      }
    }
  else
    {
    BlockAddInPlace( NL, neighborBlocks[centerBlockPos] );
    }
    //std::cout <<x<< ","<<y<<","<<z<<"\t"<< (weightsSum-neighborWeights[centerBlockPos])/(numNeighborsUsed-1);
    //std::cout << "\t" << neighborWeights[centerBlockPos] <<"\n";
    //std::cout <<x<< ","<<y<<","<<z<<"\t"<<weightsSum<<"\t"<<tmp<<"\t";
    //std::cout << numNeighborsUsed<< "\n";
}

TypedArray* 
FilterVolume::CoupeFilter
( const UniformVolume* volume, 
  const int windowRadius,
  const float beta )
{
  /*  This algorithm requires that blocks be smaller
   *  than the search neighborhood / window.
   */
  assert ( COUPE_BLOCK_RADIUS < windowRadius );

  bool BLOCKWISE_NLM = true;

  const TypedArray* inputData = volume->GetData();
  if ( ! inputData ) return NULL;

  TypedArray* filtered = TypedArray::Create( inputData->GetType(), inputData->GetDataSize() );
  const int* dims = volume->GetDims();
  const int dimX = dims[AXIS_X];
  const int dimY = dims[AXIS_Y];
  const int dimZ = dims[AXIS_Z];
  const int blockRadius = COUPE_BLOCK_RADIUS;
  const int stepSize = BLOCKWISE_NLM ? 2 : 1;
  
  std::cout << "Block radius:\t" << COUPE_BLOCK_RADIUS << "\n";
  std::cout << "Window radius:\t" << windowRadius << "\n";
  std::cout << "Beta:\t\t" << beta << std::endl;
  std::cout << "Dimensions:\t" 
            << dimX << " x " 
            << dimY << " x " 
            << dimZ << std::endl;
  
  std::vector< std::vector<Types::DataItem>* > NLsPerVoxel;

  /*  Initialize an array with a vector for each voxel,
   *  to store a set of NL estimates for each voxel.
   */
  for ( int i = 0; i < ( dimX * dimY * dimZ ); i++ )
    {
    std::vector<Types::DataItem>* tmp;
    NLsPerVoxel.push_back( tmp );
    NLsPerVoxel[i] = new std::vector<Types::DataItem>();
    } 
 
  /*  Compute the smoothing parameter 
   */
  //Types::DataItem varianceEst = ( EstimateNoiseVariance( inputData, dims ) );
  Types::DataItem stdDevEst = (Types::DataItem)3.70101;
  Types::DataItem varianceEst = stdDevEst * stdDevEst;
  std::cout << "varianceEst:\t" << varianceEst << std::endl;
  Types::DataItem smoothingParam = 2 * beta * varianceEst * COUPE_BLOCK_SIZE;
  //Types::DataItem smoothingParam = pow( 2.19037, 2 );
  std::cout << "h:\t\t"<< sqrt(smoothingParam) << std::endl;
  
  /*  Figure out where the blocks will be
   */
  TypedArray* blockLocations = TypedArray::Create(  TYPE_INT, inputData->GetDataSize() );
  int blockCount = 0;
  for ( int z = blockRadius; z < dimZ - blockRadius; z++ )
    for ( int y = blockRadius; y < dimY - blockRadius ; y++ )
      for ( int x = blockRadius; x < dimX - blockRadius; x++ )
	{
        int offset = x + dimX * ( y + dimY * z );
        if ( ( ( x % stepSize ) == 0 ) && ( ( y % stepSize ) == 0 ) && ( ( z % stepSize ) == 0 ) ) 
          {
          blockLocations->Set( 1.0, offset );
          blockCount++;
          }
        else
          {
          blockLocations->Set( 0.0, offset );
          }
        }
  std::cout << "Block count:\t" << blockCount << std::endl;


  /*  Precompute the local means and local variances maps
   */
  CoupeBlock curBlock;
  TypedArray* localMeansMap = TypedArray::Create( TYPE_DOUBLE, inputData->GetDataSize() );
  TypedArray* localVariancesMap = TypedArray::Create( TYPE_DOUBLE, inputData->GetDataSize() );
  for ( int z = blockRadius; z < dimZ - blockRadius; z++ )
    for ( int y = blockRadius; y < dimY - blockRadius ; y++ )
      for ( int x = blockRadius; x < dimX - blockRadius; x++ )
	{
        int offset = x + dimX * ( y + dimY * z );
        FilterVolume::GetCoupeBlock( curBlock, inputData, dims, x, y, z );
        Types::DataItem mean = FilterVolume::Mean( curBlock );
        Types::DataItem variance = FilterVolume::Variance( curBlock, mean );
        localMeansMap->Set( mean, offset);
        localVariancesMap->Set( variance, offset);
        }

  /*  Loop through the image, computing NL estimates
   *  for each voxel.
   */
//  CoupeBlock centerBlock;
//  CoupeBlock curNL;
  Types::DataItem blockAtCurVox = 0.0;
  for ( int Cz = blockRadius; Cz < dimZ - blockRadius; Cz ++ )
    {
    for ( int Cy = blockRadius; Cy < dimY - blockRadius; Cy ++ )
      for ( int Cx = blockRadius; Cx < dimX - blockRadius; Cx ++ ) 
	{
        int windowOffset = Cx + dimX * ( Cy + dimY * Cz );
        blockLocations->Get( blockAtCurVox, windowOffset ); 
        if ( blockAtCurVox ) // if there's a block here
          {
         
          Types::DataItem maxWeight = 0.0; // to hold the highest neighbor-weight below
          Types::DataItem weightsSum = 0.0; // to hold the sum of the neighbor weights below
          Types::DataItem outputTmp = 0.0; // to accumulate the neighbor weights in voxel-wise case
          Types::DataItem NLblock[COUPE_BLOCK_SIZE]; // to accumulate the weighted sum of neighbor blocks in block-wise case
          Types::DataItem flattenedBlock[COUPE_BLOCK_SIZE]; // to store 1-D representation of neighbor block in block-wise case
          
          /*  Initialize accumulator block for current neighborhood
           */
          if ( BLOCKWISE_NLM )
            for ( int i = 0; i < COUPE_BLOCK_SIZE; i++ ) NLblock[i] = 0.0;

          /*  Gather statistics on current voxel for use in the
           *  similarity threshold below
           */
          Types::DataItem centerMean, centerVariance;
          localMeansMap->Get( centerMean, windowOffset );
          localVariancesMap->Get( centerVariance, windowOffset );

          /*  Iterate through the blocks of the window centered at Cx,Cy,Cz
           *  ( Nx, Ny, Nz are the x,y,z coordinates of the current neighbor block )
           */
          for ( int Nz = std::max( Cz - windowRadius, 0 ); Nz < std::min( Cz + windowRadius, dimZ ); Nz++ )
            for ( int Ny = std::max( Cy - windowRadius, 0 ); Ny < std::min( Cy + windowRadius, dimZ ); Ny++ )
              for ( int Nx = std::max( Cx - windowRadius, 0 ); Nx < std::min( Cx + windowRadius, dimZ ); Nx++ )
                {
                Types::DataItem neighborDist = 0.0;  // Will hold Euclidean dist. of center and current neighbor intensities
                int neighborOffset = Nx + dimX * ( Ny + dimY * Nz );
                blockLocations->Get( blockAtCurVox, neighborOffset ); 
                if ( blockAtCurVox 
                    && ( ( Nx != Cx) || ( Ny != Cy ) || ( Nz != Cz ) ) ) // skip center block for now
                  {
                  /*  Only use blocks falling within the similarity range.
                   *  So here, compute the ratio of the means of the center
                   *  block and the neighbor block, and the ratio of their
                   *  variances.  Then proceed only if those values fall
                   *  within the desired range.
                   */
                  Types::DataItem nbMean, nbVariance;
                  localMeansMap->Get( nbMean, neighborOffset );
                  localVariancesMap->Get( nbVariance, neighborOffset );
                  Types::DataItem ratioMeans = centerMean / nbMean;
                  Types::DataItem ratioVariance = centerVariance / nbVariance;
                  Types::DataItem mu1 = (Types::DataItem)0.5;
				  Types::DataItem sigma1 = (Types::DataItem)0.95;
                  if ( ( ratioMeans > mu1 ) && ( ratioMeans <= 1 / mu1 )
                       && ( ratioVariance > sigma1 ) && ( ratioVariance <= 1 / sigma1 ) )
                    {
                   
                    /*  Initialize 1-D array for applicaton of weight to current block
                     */
                    if ( BLOCKWISE_NLM )
                      for ( int i = 0; i < COUPE_BLOCK_SIZE; i++ ) flattenedBlock[i] = 0.0;

                    /*  Iterate through center block and neighbor block in parallel,
                     *  calculating the squared Euclidean distance between the two.
                     *  Only use blocks that fall completely within the image.
                     *  ( Ix,Iy,Iz are the x,y,z coordinates of the center block voxels,
                     *    Jx,Jx,jZ are the x,y,z coordinates of the neighbor block voxels. )
                     */
                    int i = 0;
                    for ( int Bz = -blockRadius; Bz <= blockRadius; Bz++ )
                      {
                      int Iz = Cz + Bz; int Jz = Nz + Bz;
                      if ( ( Iz >= 0 ) && ( Iz < dimZ ) && ( Jz >= 0 ) && ( Jz < dimZ ) ) 
                      for ( int By = -blockRadius; By <= blockRadius; By++ )
                        {
                        int Iy = Cy + By; int Jy = Ny + By;
                        if ( ( Iy >= 0 ) && ( Iy < dimY ) && ( Jy >= 0 ) && ( Jy < dimY ) ) 
                        for ( int Bx = -blockRadius; Bx <= blockRadius; Bx++ )
                          {
                          int Ix = Cx + Bx; int Jx = Nx + Bx;
                          if ( ( Ix >= 0 ) && ( Ix < dimX ) && ( Jx >= 0 ) && ( Jx < dimX ) ) 
                            {
                            int BIoffset = Ix + dimX * ( Iy + dimY * Iz );
                            int BJoffset = Jx + dimX * ( Jy + dimY * Jz );
                            Types::DataItem BI, BJ;
                            inputData->Get( BI, BIoffset );
                            inputData->Get( BJ, BJoffset );
                            Types::DataItem distIJ = BI - BJ;
                            neighborDist += distIJ * distIJ;
                            if ( BLOCKWISE_NLM ) flattenedBlock[i++] = BJ;
                            }
                          }
                        }
                      }
                    Types::DataItem curWeight = exp( log(neighborDist) / smoothingParam );
                    weightsSum += curWeight;
                    maxWeight = ( curWeight > maxWeight ) ? curWeight : maxWeight;
                   
                    if ( BLOCKWISE_NLM )
                      {
                      /*  Weight the neighbor block and add that 
                       *  to the accumulator block
                       */
                      for ( int ii = 0; ii < COUPE_BLOCK_SIZE; ii++)
                        {
                        NLblock[ii] += curWeight * flattenedBlock[ii];
          //std::cout << curWeight << "\t" << flattenedBlock[ii] << "\n";
                        }
                      }
                      else  // voxel-wise case
                        {
                        /*  Weight the neighbor voxel and add that 
                         *  to the accumulator value
                         */
                        Types::DataItem neighbValue;
                        inputData->Get( neighbValue, neighborOffset );
                        outputTmp += curWeight * neighbValue;
                        }

                    } // end similarity-threshold test
                  } // end test for J != I
                } // end loop through blocks of current window
          
          /*  Current center block gets weighted by
           *  the maximum neighbor-weight
           */ 
          Types::DataItem centerValue;
          inputData->Get( centerValue, windowOffset );     
          outputTmp += maxWeight * centerValue;
          weightsSum += maxWeight;

          if ( BLOCKWISE_NLM )
            {
            /*  Push the accumulator block into place for averaging  
             *  later 
             *  ( Ix,Iy,Iz are the coordinates the center block voxels )
             */
            int i = 0;
            for ( int Iz = Cz - blockRadius; Iz <= Cz + blockRadius; Iz++ )
              {
              if ( ( Iz >= 0 ) && ( Iz < dimZ ) )
              for ( int Iy = Cy - blockRadius; Iy <= Cy + blockRadius; Iy++ )
                {
                if ( ( Iy >= 0 ) && ( Iy < dimY ) )
                for ( int Ix = Cx - blockRadius; Ix <= Cx + blockRadius; Ix++ )
                  {
                  if ( ( Ix >= 0 ) && ( Ix < dimX ) )
                    {
                    int Voffset = Ix + dimX * ( Iy + dimY * Iz );
          //std::cout << Voffset << "\t" << NLblock[i] << "\n";
                    NLsPerVoxel[Voffset]->push_back( NLblock[i] );
                    i++;
                    }
                  }
                }
              }
            }
          else
            {
            /*  Set the output voxel
             */
            filtered->Set( outputTmp / weightsSum, windowOffset );
            }

          } // end if ( blockAtCurVox )
        } // end loop through pixels of image 

//                    FilterVolume::ComputeNLWithinWindow( curNL, blockLocations, inputData, dims, smoothingParam,
//                                                            Cx, Cy, Cz, 
//                                                            windowRadius, 
//                                                            beta,
//                                                            localMeansMap, 
//                                                            localVariancesMap, 
//                                                            centerBlock );
//                    /*  Push the values from curNL into the vector
//                    *  containing the NL estimates for the respective
//                    *  corresponding voxels.  (curNL contains an NL
//                    *  estimate for each voxel in this block).
//                    */
//                    int curNLIndex = 0;
//                    for ( int k = z - blockRadius; k <= z + blockRadius; ++k )
//                      if ( ( k >= 0 ) && ( k < dimZ ) )
//                        {
//                        for ( int j = y - blockRadius; j <= y + blockRadius; ++j )
//                          if ( ( j >= 0 ) && ( j < dimY ) )
//                            {
//                            for ( int i = x - blockRadius; i <= x + blockRadius; ++i )
//                              if ( ( i >= 0 ) && ( i < dimX ) )
//                                {
//                                int offset = i + dimX * ( j + dimY * k );
//                                NLsPerVoxel[offset]->push_back( curNL[curNLIndex] ); 
//                                curNLIndex++;
//                                }
//                            }
//                        }
    }



  if ( BLOCKWISE_NLM )
    {
    /*  Take the average of the NL estimates for each voxel
    */
    int NLsPerVoxelHist[40] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    for ( int z = 0; z < dimZ; ++z )
      for ( int y = 0; y < dimY; ++y )
        for ( int x = 0; x < dimX; ++x ) 
          {
          int offset = x + dimX * ( y + dimY * z );
          NLsPerVoxelHist[NLsPerVoxel.at( offset )->size()]++;
          Types::DataItem curNLSum = 0.0;
          for ( std::vector<Types::DataItem>::iterator it = NLsPerVoxel.at( offset )->begin(); it != NLsPerVoxel.at( offset )->end(); it++ )
            curNLSum += *it;

          Types::DataItem restoredVoxel = curNLSum / NLsPerVoxel.at( offset )->size();
          //std::cout << offset << "\t" << curNLSum << "\t" << NLsPerVoxel.at( offset )->size() << "\n";

          filtered->Set( restoredVoxel , offset );

          if ( ( x == 91 ) && ( y == 1 ) && ( z == 93 ) ) 
            {
            std::cout << std::endl;
            std::cout << std::endl;
            for ( std::vector<Types::DataItem>::iterator it = NLsPerVoxel.at( offset )->begin(); it != NLsPerVoxel.at( offset )->end(); it++ )
              {
              std::cout << *it << "..." << "\n";
              }
            std::cout << restoredVoxel << "." << "\n";
            }

          }

    /*  Print out the histogram of NL counts
     */ 
    std::cout << std::endl;
    std::cout << "NL-count histogram:" << std::endl;
    std::cout << "0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19" << std::endl;
    for ( int i = 0; i < 20; i++ )
      std::cout << NLsPerVoxelHist[i] << "\t";
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "20\t21\t22\t23\t24\t25\t26\t27\t28\t29\t30\t31\t32\t33\t34\t35\t36\t37\t38\t39" << std::endl;
    for ( int i = 20; i < 40; i++ )
      std::cout << NLsPerVoxelHist[i] << "\t";
    std::cout << std::endl;
    std::cout << std::endl;
    }
 

  return filtered;
}

} // namespace cmtk
