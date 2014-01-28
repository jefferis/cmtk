/*
//
//  Copyright 2010 Torsten Rohlfing
//
//  Copyright 2011, 2014 SRI International
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
//  $Revision: 3097 $
//
//  $LastChangedDate: 2011-04-06 13:07:22 -0700 (Wed, 06 Apr 2011) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include "cmtkImageOperationMatchIntensities.h"

#include <IO/cmtkVolumeIO.h>

#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkTypedArrayFunctionHistogramMatching.h>

cmtk::ImageOperationMatchIntensities::ImageOperationMatchIntensities( const Self::Mode mode, const std::string& referenceImagePath )
  : m_Mode( mode )
{
  UniformVolume::SmartPtr referenceImage = VolumeIO::Read( referenceImagePath );
  if ( ! referenceImage )
    {
    StdErr << "ERROR: cannot read image " << referenceImagePath << "\n";
    throw ExitException( 1 );
    }

  this->m_ReferenceData = referenceImage->GetData();
  if ( ! this->m_ReferenceData )
    {
    StdErr << "ERROR: read geometry but could not read pixel data from " << referenceImagePath << "\n";
    throw ExitException( 1 );
    }
}

cmtk::UniformVolume::SmartPtr  
cmtk::ImageOperationMatchIntensities::Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    TypedArray& volumeData = *(volume->GetData());
    switch ( this->m_Mode ) 
      {
      case Self::MATCH_HISTOGRAMS:
	volumeData.ApplyFunctionObject( cmtk::TypedArrayFunctionHistogramMatching( volumeData, *(this->m_ReferenceData) ) );
	break;
      case Self::MATCH_MEAN_SDEV:
      {
      cmtk::Types::DataItem rMean, rVar;
      this->m_ReferenceData->GetStatistics( rMean, rVar );
      
      cmtk::Types::DataItem mMean, mVar;
      volumeData.GetStatistics( mMean, mVar );
      
      const cmtk::Types::DataItem scale = sqrt( rVar / mVar );
      const cmtk::Types::DataItem offset = rMean - scale * mMean;
      volumeData.Rescale( scale, offset );
      break;
      }
      }
    return volume;
  }
  
void
cmtk::ImageOperationMatchIntensities::NewMatchHistograms( const char* referenceImagePath )
{
  ImageOperation::m_ImageOperationList.push_back( SmartPtr( new Self( Self::MATCH_HISTOGRAMS, referenceImagePath ) ) );
}
  
void
cmtk::ImageOperationMatchIntensities::NewMatchMeanSDev( const char* referenceImagePath )
{
  ImageOperation::m_ImageOperationList.push_back( SmartPtr( new Self( Self::MATCH_MEAN_SDEV, referenceImagePath ) ) );
}

