/*
//
//  Copyright 2003 Calvin R. Maurer, Jr.
//
//  Copyright 1997-2009 Torsten Rohlfing
//
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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkUniformDistanceMap.h>

#include <cmtkDataTypeTraits.h>

#include <cmtkThreadPool.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class TDistanceDataType>
UniformDistanceMap<TDistanceDataType>
::UniformDistanceMap
( const UniformVolume& volume, const byte flags, const Types::DataItem value, const Types::DataItem window )
  : Superclass( volume.m_Dims, volume.Size )
{
  this->BuildDistanceMap( volume, flags, value, window );

  this->m_IndexToPhysicalMatrix = volume.m_IndexToPhysicalMatrix;

  this->SetOffset( volume.m_Offset );
  this->m_MetaInformation = volume.m_MetaInformation;
}

template<class TDistanceDataType>
void
UniformDistanceMap<TDistanceDataType>
::BuildDistanceMap
( const UniformVolume& volume, const byte flags, const Types::DataItem value, const Types::DataItem window )
{
  TypedArray::SmartPtr distanceArray = TypedArray::SmartPtr( TypedArray::Create( DataTypeTraits<DistanceDataType>::DataTypeID, volume.GetNumberOfPixels() ) );
  DistanceDataType *Distance = static_cast<DistanceDataType*>( distanceArray->GetDataPtr() );

  byte inside = ( flags & UniformDistanceMap::INSIDE ) ? 0 : 1;
  byte outside = 1 - inside;

  const TypedArray* Feature = volume.GetData();

  Types::DataItem c;
  DistanceDataType *p = Distance;
  if ( flags & UniformDistanceMap::VALUE_EXACT ) 
    {
    for ( size_t i = 0; i < volume.GetNumberOfPixels(); i++, p++ ) 
      {
      if ( Feature->Get( c, i ) )
	{
	*p = (c == value) ? inside : outside;
	}
      else
	{
	*p = outside;
	}
      }
    } 
  else if ( flags & UniformDistanceMap::VALUE_THRESHOLD ) 
    {
    for ( size_t i = 0; i < volume.GetNumberOfPixels(); i++, p++ ) 
      {
      if ( Feature->Get( c, i ) )
	{
	*p = (c >= value) ? inside : outside;
	}
      else
	{
	*p = outside;
	}
      }
    } 
  else if ( flags & UniformDistanceMap::VALUE_WINDOW ) 
    {
    for ( size_t i = 0; i < volume.GetNumberOfPixels(); i++, p++ ) 
      {
      if ( Feature->Get( c, i ) )
	{
	*p = (fabs(c - value)<=window) ? inside : outside;
	}
      else
	{
	*p = outside;
	}
      }
    } 
  else
    {
    for ( size_t i = 0; i < volume.GetNumberOfPixels(); i++, p++ ) 
      {
      if ( Feature->Get( c, i ) )
	{
	*p = (c) ? inside : outside;
	}
      else
	{
	*p = outside;
	}
      }
    }
  
  this->ComputeEDT( Distance );

  p = Distance;
  for ( size_t i = 0; i < volume.GetNumberOfPixels(); ++i, ++p ) 
    {
#if defined(_MSC_VER) || defined(__SUNPRO_CC)
    *p = static_cast<DistanceDataType>( sqrt( (double)*p ) );
#else
    *p = static_cast<DistanceDataType>( sqrt( *p ) );
#endif
    }
  this->SetData( distanceArray );
}

template<class TDistanceDataType>
void
UniformDistanceMap<TDistanceDataType>
::ComputeEDT( DistanceDataType *const distance )
{
  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  const size_t numberOfThreads = threadPool.GetNumberOfThreads();
  const size_t numberOfTasks = 4 * numberOfThreads - 3;
  
  this->m_G.resize( numberOfThreads );
  this->m_H.resize( numberOfThreads );
  
  std::vector<typename Self::ThreadParametersEDT> params( numberOfTasks );
  for ( size_t idx = 0; idx < numberOfTasks; ++idx )
    {
    params[idx].thisObject = this;
    params[idx].m_Distance = distance;
    }

  threadPool.Run( ComputeEDTThreadPhase1, params );
  threadPool.Run( ComputeEDTThreadPhase2, params );
}

template<class TDistanceDataType>
void
UniformDistanceMap<TDistanceDataType>
::ComputeEDTThreadPhase1
( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  ThreadParametersEDT* params = static_cast<ThreadParametersEDT*>( args );
  Self* This = params->thisObject;
  const Self* ThisConst = This;

  /* nXY is number of voxels in each plane (xy) */
  /* nXYZ is number of voxels in 3D image */
  const size_t nXY = ThisConst->m_Dims[0] * ThisConst->m_Dims[1];
  
  /* compute D_2 */
  /* call edtComputeEDT_2D for each plane */
  DistanceDataType *p = params->m_Distance + nXY * taskIdx;
  for ( int k = taskIdx; k < ThisConst->m_Dims[2]; k += taskCnt, p += nXY * taskCnt ) 
    {
    This->ComputeEDT2D( p, This->m_G[threadIdx], This->m_H[threadIdx] );
    }
}

template<class TDistanceDataType>
void
UniformDistanceMap<TDistanceDataType>
::ComputeEDTThreadPhase2
( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t )
{
  ThreadParametersEDT* params = static_cast<ThreadParametersEDT*>( args );
  Self* This = params->thisObject;
  const Self* ThisConst = This;

  const size_t nXY = ThisConst->m_Dims[0] * ThisConst->m_Dims[1];
  /* compute D_3 */
  /* solve 1D problem for each column (z direction) */
  std::vector<DistanceDataType> f( This->m_Dims[2] );

  for ( size_t i = taskIdx; i < nXY; i += taskCnt ) 
    {
    /* fill array f with D_2 distances in column */
    /* this is essentially line 4 in Procedure VoronoiEDT() in tPAMI paper */
    DistanceDataType *p = params->m_Distance + i;
    DistanceDataType *q = &f[0];
    for ( int k = 0; k < ThisConst->m_Dims[2]; k++, p += nXY, q++) 
      {
      *q = *p;
      }
    
    /* call edtVoronoiEDT */
    if ( This->VoronoiEDT( &f[0], ThisConst->m_Dims[2], static_cast<DistanceDataType>( ThisConst->m_Delta[2] ), This->m_G[threadIdx], This->m_H[threadIdx] ) ) 
      {
      p = params->m_Distance + i;
      DistanceDataType *q = &f[0];
      for ( int k = 0; k < ThisConst->m_Dims[2]; k++, p += nXY, q++ ) 
	{
	*p = *q;
	}
      }
    }
}

template<class TDistanceDataType>
void
UniformDistanceMap<TDistanceDataType>
::ComputeEDT2D
( DistanceDataType *const plane, std::vector<DistanceDataType>& gTemp, std::vector<DistanceDataType>& hTemp )
  /*
   * This procedure computes the squared EDT of a 2D binary image with
   * anisotropic voxels. See notes for edtComputeEDT_2D. The difference 
   * relative to edtComputeEDT_2D is that the edt is a float array instead of 
   * a long array, and there are additional parameters for the image voxel 
   * dimensions wX and wY.
   */
{  
  /* compute D_1 as simple forward-and-reverse distance propagation */
  /* (instead of calling edtVoronoiEDT) */
  /* D_1 is distance to closest feature voxel in row (x direction) */
  /* it is possible to use a simple distance propagation for D_1  because */
  /* L_1 and L_2 norms are equivalent for 1D case */
  DistanceDataType *p;
  for ( int j = 0; j < this->m_Dims[1]; j++ ) 
    {
    /* forward pass */
    p = plane + j * this->m_Dims[0];
    DistanceDataType d = static_cast<DistanceDataType>( EDT_MAX_DISTANCE_SQUARED );
    for ( int i = 0; i < this->m_Dims[0]; i++, p++ ) 
      {
      /* set d = 0 when we encounter a feature voxel */
      if ( *p ) 
	{
	*p = d = 0;
	}
      /* increment distance ... */
      else 
	if ( d != EDT_MAX_DISTANCE_SQUARED ) 
	  {
	  *p = ++d;
	  }
      /* ... unless we haven't encountered a feature voxel yet */
	else
	  {
	  *p = static_cast<DistanceDataType>( EDT_MAX_DISTANCE_SQUARED );
	  }
      }
    
    /* reverse pass */
    if ( *(--p) != EDT_MAX_DISTANCE_SQUARED ) 
      {
      DistanceDataType d = static_cast<DistanceDataType>( EDT_MAX_DISTANCE_SQUARED );
      for ( int i = this->m_Dims[0] - 1; i >= 0; i--, p-- ) 
	{
	/* set d = 0 when we encounter a feature voxel */
	if ( *p == 0 ) 
	  {
	  d = 0;
	  }
	/* increment distance after encountering a feature voxel */
	else
	  if ( d != EDT_MAX_DISTANCE_SQUARED ) 
	    {
	    /* compare forward and reverse distances */
	    if ( ++d < *p ) 
	      {
	      *p = d;
	      }
	    }
	
	/* square distance */
	/* (we use squared distance in rest of algorithm) */
	*p = static_cast<DistanceDataType>( *p * m_Delta[0] );
	*p *= *p;
	}
      }
    }
  
  /* compute D_2 = squared EDT */
  /* solve 1D problem for each column (y direction) */
  std::vector<DistanceDataType> f( this->m_Dims[1] );
  for ( int i = 0; i < this->m_Dims[0]; i++ ) 
    {
    /* fill array f with D_1 distances in column */
    /* this is essentially line 4 in Procedure VoronoiEDT() in tPAMI paper */
    DistanceDataType *p = plane + i;
    DistanceDataType *q = &f[0];
    for ( int j = 0; j < this->m_Dims[1]; j++, p += this->m_Dims[0], q++) 
      {
      *q = *p;
      }
    
    /* call edtVoronoiEDT */
    if ( this->VoronoiEDT( &f[0], this->m_Dims[1], static_cast<DistanceDataType>( this->m_Delta[1] ), gTemp, hTemp  ) ) 
      {
      DistanceDataType *p = plane + i;
      DistanceDataType *q = &f[0];
      for ( int j = 0; j < this->m_Dims[1]; j++, p += this->m_Dims[0], q++ ) 
	{
	*p = *q;
	}
      }
    }
  
} /* edtComputeEDT_2D_anisotropic */

template<class TDistanceDataType>
bool
UniformDistanceMap<TDistanceDataType>
::VoronoiEDT
( DistanceDataType *const distanceSoFar, const int nSize, const DistanceDataType delta,
  std::vector<DistanceDataType>& gTemp, std::vector<DistanceDataType>& hTemp )
{
  long i, l, n_S;
  DistanceDataType a, b, c, v, lhs, rhs;
  
  /* alloc arrays if this is first call to procedure, or if arrays */
  /* are too small and need to be reallocated */
  gTemp.resize( nSize );
  hTemp.resize( nSize );
  
  DistanceDataType* g = &(gTemp[0]);
  DistanceDataType* h = &(hTemp[0]);

  /* construct partial Voronoi diagram */
  /* this loop is lines 1-14 in Procedure edtVoronoiEDT() in tPAMI paper */
  /* note we use 0 indexing in this program whereas paper uses 1 indexing */
  DistanceDataType deltai = 0;
  for (i = 0, l = -1; i < nSize; i++, deltai += delta) 
    {
    /* line 4 */
    if ( distanceSoFar[i] != EDT_MAX_DISTANCE_SQUARED ) 
      {
      /* line 5 */
      if ( l < 1 ) 
	{
	/* line 6 */
	g[++l] = distanceSoFar[i];
	h[l] = deltai;
	}
      /* line 7 */
      else 
	{
	/* line 8 */
	while (l >= 1) 
	  {
	  /* compute removeEDT() in line 8 */
	  v = h[l];
	  a = v - h[l-1];
	  b = deltai - v;
	  c = a + b;
	  /* compute Eq. 2 */
	  if ((c*g[l] - b*g[l-1] - a*distanceSoFar[i] - a*b*c) > 0) 
	    {
	    /* line 9 */
	    l--;
	    } 
	  else 
	    {
	    break;
	    }
	  }
	/* line 11 */
	g[++l] = distanceSoFar[i];
	h[l] = deltai;
	}
      }
    }
  /* query partial Voronoi diagram */
  /* this is lines 15-25 in Procedure edtVoronoiEDT() in tPAMI paper */
  /* lines 15-17 */
  if ((n_S = l + 1) == 0) 
    {
    return false;
    }
  
  /* lines 18-19 */
  deltai = 0;
  for (i = 0, l = 0; i < nSize; i++, deltai += delta) 
    {
    /* line 20 */
    /* we reduce number of arithmetic operations by taking advantage of */
    /* similarities in successive computations instead of treating them as */
    /* independent ones */
    a = h[l] -  deltai;
    lhs = g[l] + a * a;
    while (l < n_S - 1)
      {
      a = h[l+1] - deltai;
      rhs = g[l+1] + a * a;
      if (lhs > rhs) 
	{
	/* line 21 */
	l++;
	lhs = rhs;
	} 
      else
	{
	break;
	}
      }
    
    /* line 23 */
    /* we put distance into the 1D array that was passed; */
    /* must copy into EDT in calling procedure */
    distanceSoFar[i] = lhs;
    }
  
  /* line 25 */
  /* return 1 if we queried diagram, 0 if we returned because n_S = 0 */
  return true;
}

//@}

} // namespace cmtk

template class cmtk::UniformDistanceMap<float>;
template class cmtk::UniformDistanceMap<double>;
template class cmtk::UniformDistanceMap<long int>;
