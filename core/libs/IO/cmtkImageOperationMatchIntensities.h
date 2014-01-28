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

#ifndef __cmtkImageOperationMatchIntensities_h_included_
#define __cmtkImageOperationMatchIntensities_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: match intensities to another image.
class ImageOperationMatchIntensities
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// This class.
  typedef ImageOperationMatchIntensities Self;

  /// Superclass.
  typedef ImageOperation Superclass;

  /// Operation mode.
  typedef enum
  {
    /// Match histograms.
    MATCH_HISTOGRAMS,
    /// Match mean and standard deviation
    MATCH_MEAN_SDEV
  } Mode;

  /// Constructor.
  ImageOperationMatchIntensities( const Self::Mode mode /*!< Operation mode.*/, const std::string& referenceImagePath /*!< Path of the reference image to match intensites to.*/ );

  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume );
  
  /// Create new operation to match histograms.
  static void NewMatchHistograms( const char* referenceImagePath /*!< Path of the reference image to match intensites to.*/ );
  
  /// Create new operation to match image mean and standard deviation.
  static void NewMatchMeanSDev( const char* referenceImagePath /*!< Path of the reference image to match intensites to.*/ );
  
private:
  /// Operation mode.
  Self::Mode m_Mode;

  /// Reference image data.
  TypedArray::SmartPtr m_ReferenceData;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationMatchIntensities_h_included_
