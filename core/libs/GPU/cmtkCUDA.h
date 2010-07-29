/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkCUDA_h_included_
#define __cmtkCUDA_h_included_

#include <cmtkconfig.h>

#include <cuda_runtime_api.h>

#include <cstdio>
#include <cstdlib>

#define cmtkCheckCallCUDA(cmd) \
  { const cudaError_t cudaError = cmd; if ( cudaError != cudaSuccess ) { fprintf( stderr, "CUDA command failed with error '%s' at %s:%d\n", cudaGetErrorString( cudaError ), __FILE__, __LINE__ ); exit(1); } }
  
#define cmtkCheckLastErrorCUDA \
  { const cudaError_t cudaError = cudaGetLastError(); if ( cudaError != cudaSuccess ) { fprintf( stderr, "CUDA error '%s' at %s:%d\n", cudaGetErrorString( cudaError ), __FILE__, __LINE__ ); exit( 1 ); } }


#endif // #ifndef __cmtkCUDA_h_included_
