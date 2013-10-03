/*
//
//  Copyright 2010-2013 SRI International
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

#include "cmtkSimpleLevelset.h"

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkUniformVolumePainter.h>
#include <Base/cmtkUniformVolumeGaussianFilter.h>
#include <Base/cmtkUnits.h>

#include <System/cmtkProgress.h>
#include <System/cmtkDebugOutput.h>

#ifdef CMTK_BUILD_DEMO
#  include <IO/cmtkVolumeIO.h>
#endif // #ifdef CMTK_BUILD_DEMO

void
cmtk::SimpleLevelset
::InitializeCenteredSphere()
{
  this->m_Levelset = this->m_Volume->CloneGrid();
  this->m_Levelset->CreateDataArray( TYPE_FLOAT );
  this->m_Levelset->GetData()->Fill( -1.0 );
  
  FixedVector<3,int> center( this->m_Volume->GetDims() );
  center /= 2;
  
  UniformVolumePainter painter( this->m_Levelset );
  painter.DrawSphere( center, this->m_ScaleInitialSphere * ((this->m_Levelset->GetDims()[0]+this->m_Levelset->GetDims()[1]+this->m_Levelset->GetDims()[2])/6), 1.0 );
}

void
cmtk::SimpleLevelset
::Evolve( const int numberOfIterations, const bool forceIterations )
{
  const size_t numberOfPixels = this->m_Volume->GetNumberOfPixels();
  size_t nInsideOld = 0, nInside = 1;

  Progress::Begin( 0, numberOfIterations, 1, "Levelset Evolution" );
  for ( int it = 0; (it < numberOfIterations) && ((nInside!=nInsideOld) || forceIterations); ++it )
    {
    Progress::SetProgress( it );

    nInsideOld = nInside;
    nInside = 0;
    Types::DataItem insideSum = 0, outsideSum = 0;

    this->m_Levelset->SetData( UniformVolumeGaussianFilter( this->m_Levelset ).GetFiltered3D( this->m_FilterSigma ) );
#pragma omp parallel for reduction(+:nInside) reduction(+:insideSum) reduction(+:outsideSum)
    for ( int n = 0; n < static_cast<int>( numberOfPixels ); ++n )
      {
      if ( this->m_Levelset->GetDataAt( n ) > 0 )
	{
	insideSum += this->m_Volume->GetDataAt( n );
	++nInside;
	}
      else
	outsideSum += this->m_Volume->GetDataAt( n );
      }

    if ( nInside == 0 )
      throw Self::DegenerateLevelsetException();

    const size_t nOutside = numberOfPixels - nInside;

    if ( nOutside == 0 )
      throw Self::DegenerateLevelsetException();

    const Types::DataItem ratioInOut = 1.0 * nInside / nOutside;
    
    const Types::DataItem mInside = insideSum / nInside;
    const Types::DataItem mOutside = outsideSum / nOutside;

    DebugOutput( 1 ) << it << " IN: " << nInside << "  " << mInside << "  OUT: " << nOutside << "  " << mOutside << "\r";
    
#pragma omp parallel for
    for ( int n = 0; n < static_cast<int>( numberOfPixels ); ++n )
      {
      const Types::DataItem data = this->m_Volume->GetDataAt( n );
      const Types::DataItem zInside = fabs( mInside - data );
      const Types::DataItem zOutside = fabs( mOutside - data );
      Types::DataItem newLevel = this->m_Levelset->GetDataAt( n );
      if ( zInside>zOutside )
	{
	newLevel -= this->m_TimeDelta * ratioInOut;
	}
      else
	{
	newLevel += this->m_TimeDelta / ratioInOut;
	}
      this->m_Levelset->SetDataAt( std::min<Types::DataItem>( this->m_LevelsetThreshold, std::max<Types::DataItem>( -this->m_LevelsetThreshold, newLevel ) ), n );
      }

#ifdef CMTK_BUILD_DEMO
    UniformVolume::SmartConstPtr slice = this->m_Levelset->ExtractSlice( AXIS_Z, this->m_Levelset->m_Dims[2] / 2 );
    
    char path[PATH_MAX];
    snprintf( path, PATH_MAX, "levelset-%03d.nii", it );
    VolumeIO::Write( *slice, path );
#endif
    }

  Progress::Done();
}

cmtk::UniformVolume::SmartPtr&
cmtk::SimpleLevelset
::GetLevelset( const bool binarize, const float threshold )
{
  if ( binarize )
    {
    this->m_Levelset->GetData()->Binarize( threshold );
    this->m_Levelset->SetData( TypedArray::SmartPtr( this->m_Levelset->GetData()->Convert( TYPE_BYTE ) ) );
    }

  return this->m_Levelset;
}
