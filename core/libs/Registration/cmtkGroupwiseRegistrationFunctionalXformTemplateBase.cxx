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

//#define DEBUG_COMM

#include "Registration/cmtkGroupwiseRegistrationFunctionalXformTemplateBase.h"

#include "Base/cmtkMathUtil.h"
#include "IO/cmtkVolumeIO.h"

#ifdef CMTK_BUILD_MPI
#  include <mpi.h>
#  include "IO/cmtkMPI.h"
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TXform>
GroupwiseRegistrationFunctionalXformTemplateBase<TXform>::GroupwiseRegistrationFunctionalXformTemplateBase() :
  m_HistogramBins( 64 ),
  m_HistogramKernelRadiusMax( 0 ),
  m_MaxRelativeNumberOutsidePixels( 0.99 ), // if there is an image with more then 99% pixels outside FOV, registration probably failed
  m_CropImageHistograms( false )
{}

template<class TXform>
GroupwiseRegistrationFunctionalXformTemplateBase<TXform>::~GroupwiseRegistrationFunctionalXformTemplateBase()
{
}

template<class TXform>
void
GroupwiseRegistrationFunctionalXformTemplateBase<TXform>
::SetNumberOfHistogramBins( const size_t numberOfHistogramBins )
{
  this->m_HistogramBins = numberOfHistogramBins;
  if ( this->m_OriginalImageVector.size() )
    {
    std::cerr << "WARNING: you called GroupwiseRegistrationFunctionalBase::SetNumberOfHistogramBins(),\n"
	      << "         but target images were already set. To be safe, I am re-generating\n"
	      << "         pre-scaled images.\n\n";
    this->SetTargetImages( this->m_OriginalImageVector );
    }
}

template<class TXform>
UniformVolume*
GroupwiseRegistrationFunctionalXformTemplateBase<TXform>
::PrepareSingleImage( UniformVolume::SmartPtr& image )
{
  UniformVolume* newTargetImage = this->Superclass::PrepareSingleImage( image );

  TypedArray::SmartPtr data = newTargetImage->GetData();
  if ( this->m_CropImageHistograms )
    {
    data->PruneHistogram( true, false, this->m_HistogramBins );
    }
  
  data->Rescale( 1.0 * (this->m_HistogramBins-1) / data->GetRange().Width(), this->m_HistogramKernelRadiusMax );
  
  newTargetImage->SetData( TypedArray::SmartPtr( data->Convert( TYPE_BYTE ) ) );
  return newTargetImage;
}

template<class TXform>
void
GroupwiseRegistrationFunctionalXformTemplateBase<TXform>
::PrepareTargetImages()
{
  this->m_ImageVector.resize( this->m_OriginalImageVector.size() );

#ifdef CMTK_BUILD_MPI
  // using MPI, prepare only some of the images locally, obtain others from other nodes
  const size_t imageFrom = this->m_RankMPI;
  const size_t imageSkip = this->m_SizeMPI;
#else
  const size_t imageFrom = 0;
  const size_t imageSkip = 1;
#endif
  for ( size_t i = imageFrom; i < this->m_ImageVector.size(); i += imageSkip )
    {
    this->m_ImageVector[i] = UniformVolume::SmartPtr( this->PrepareSingleImage( this->m_OriginalImageVector[i] ) );
    }
  
#ifdef CMTK_BUILD_MPI
  // obtain filtered, scaled image data from other nodes
  for ( size_t i = 0; i < this->m_ImageVector.size(); ++i )
    {
    cmtk::mpi::Broadcast( MPI::COMM_WORLD, this->m_ImageVector[i], i % this->m_SizeMPI );
    }
#endif

  this->m_PrivateUserBackgroundValue = this->m_UserBackgroundValue + this->m_HistogramKernelRadiusMax;
}

//@}

} // namespace cmtk

#include "Base/cmtkAffineXform.h"
#include "Base/cmtkSplineWarpXform.h"

template class cmtk::GroupwiseRegistrationFunctionalXformTemplateBase<cmtk::AffineXform>;
template class cmtk::GroupwiseRegistrationFunctionalXformTemplateBase<cmtk::SplineWarpXform>;
