/*
//
//  Copyright 2009-2012 SRI International
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

#include "cmtkImageOperationOtsuThreshold.h"

#include <Base/cmtkHistogramOtsuThreshold.h>

#include <System/cmtkDebugOutput.h>

cmtk::UniformVolume::SmartPtr
cmtk::ImageOperationOtsuThreshold::Apply( cmtk::UniformVolume::SmartPtr& volume )
{
  cmtk::TypedArray& volumeData = *(volume->GetData());

  const Types::DataItem threshold = HistogramOtsuThreshold< Histogram<unsigned int> >( *(volumeData.GetHistogram( this->m_NumberOfBins )) ).Get();

  DebugOutput( 5 ) << "INFO: Otsu binarization threshold = " << threshold << "\n";

  volumeData.Binarize( threshold );

  return volume;
}
