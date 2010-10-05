/*
//
//  Copyright 2009-2010 SRI International
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

#ifndef __cmtkImageOperationGaussFilter_h_included_
#define __cmtkImageOperationGaussFilter_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>
#include <Base/cmtkUniformVolumeFilter.h>

namespace
cmtk
{

/// Image operation: grid downsampling.
class ImageOperationGaussFilter
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationGaussFilter( const Units::GaussianSigma& sigma ) : m_Sigma( sigma ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    volume->SetData( UniformVolumeFilter( volume ).GetDataGaussFiltered( this->m_Sigma, 0.001 /* kernel truncation approximation error threshold */ ) );
    return volume;
  }
  
  /// Create a new filter based on sigma parameter.
  static void NewSigma( const double sigma )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationGaussFilter( Units::GaussianSigma( sigma ) ) ) );
  }
  
  /// Create a new filter based on full-width-at-half-maximum parameter.
  static void NewFWHM( const double fwhm )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationGaussFilter( Units::GaussianFWHM( fwhm ) ) ) );
  }
  
private:
  /// Kernel with specified by coefficient sigma.
  Units::GaussianSigma m_Sigma;
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationGaussFilter_h_included_
