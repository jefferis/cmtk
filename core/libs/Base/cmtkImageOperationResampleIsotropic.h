/*
//
//  Copyright 2009-2011, 2013-2014 SRI International
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

#ifndef __cmtkImageOperationResampleIsotropic_h_included_
#define __cmtkImageOperationResampleIsotropic_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: anisotropic resampling
class ImageOperationResampleIsotropic
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// This class.
  typedef ImageOperationResampleIsotropic Self;

  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr Apply( cmtk::UniformVolume::SmartPtr& volume );
  
  /// Create a new resampler.
  static void New( const double resolution );
  
  /// Create a new resampler.
  static void NewExact( const double resolution );
  
private:
  /// Constructor:
  ImageOperationResampleIsotropic( const double resolution, const bool exact = false ) : m_Resolution( resolution ), m_Exact( exact ) {}
  
  /// Anisotropic resampling resolution
  double m_Resolution;

  /// Flag for exact vs. approximate resampling.
  bool m_Exact;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationResampleIsotropic_h_included_
