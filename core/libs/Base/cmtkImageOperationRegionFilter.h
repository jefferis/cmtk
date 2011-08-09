/*
//
//  Copyright 2009-2011 SRI International
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

#ifndef __cmtkImageOperationRegionFilter_h_included_
#define __cmtkImageOperationRegionFilter_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>
#include <Base/cmtkDataGridFilter.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: region filtering.
class ImageOperationRegionFilter
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// This class.
  typedef ImageOperationRegionFilter Self;

  /// Types of filters supported by this class.
  typedef enum
  {
    MEDIAN, 
    MEAN,
    FAST_MEAN,
    VARIANCE,
    FAST_VARIANCE,
    THIRD_MOMENT,
    STANDARD_DEVIATION,
    SMOOTHNESS
  } OperatorEnum;

  /// Constructor:
  ImageOperationRegionFilter( const Self::OperatorEnum op, const int radiusX, const int radiusY, const int radiusZ ) 
    : m_Operator( op ), m_RadiusX( radiusX ), m_RadiusY( radiusY ), m_RadiusZ( radiusZ ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume );
  
  /// Create a new region median filter operation.
  static void NewMedian( const char* arg /*!< Region size argument: either "XYZ" or "X,Y,Z" */ )
  {
    Self::NewGeneric( Self::MEDIAN, arg );
  }
  
  /// Create a new region mean filter operation.
  static void NewMean( const char* arg /*!< Region size argument: either "XYZ" or "X,Y,Z" */ )
  {
    Self::NewGeneric( Self::MEAN, arg );
  }
  
  /// Create a new fast region mean filter operation.
  static void NewFastMean( const char* arg /*!< Region size argument: either "XYZ" or "X,Y,Z" */ )
  {
    Self::NewGeneric( Self::FAST_MEAN, arg );
  }
  
  /// Create a new region variance filter operation.
  static void NewVariance( const char* arg /*!< Region size argument: either "XYZ" or "X,Y,Z" */ )
  {
    Self::NewGeneric( Self::VARIANCE, arg );
  }
  
  /// Create a new fast region variance filter operation.
  static void NewFastVariance( const char* arg /*!< Region size argument: either "XYZ" or "X,Y,Z" */ )
  {
    Self::NewGeneric( Self::FAST_VARIANCE, arg );
  }
  
  /// Create a new region third moment filter operation.
  static void NewThirdMoment( const char* arg /*!< Region size argument: either "XYZ" or "X,Y,Z" */ )
  {
    Self::NewGeneric( Self::THIRD_MOMENT, arg );
  }
  
  /// Create a new region standard deviation filter operation.
  static void NewStandardDeviation( const char* arg /*!< Region size argument: either "XYZ" or "X,Y,Z" */ )
  {
    Self::NewGeneric( Self::STANDARD_DEVIATION, arg );
  }
  
  /// Create a new region smoothness filter operation.
  static void NewSmoothness( const char* arg /*!< Region size argument: either "XYZ" or "X,Y,Z" */ )
  {
    Self::NewGeneric( Self::SMOOTHNESS, arg );
  }
  
private:
  /// Parse region size argument and create new filter object.
  static void NewGeneric( const Self::OperatorEnum op /*!< The operation to perform by the new object.*/, const char* arg /*!< Region size argument: either "XYZ" or "X,Y,Z" */ );
  
  /// The operator this object will apply to its input.
  Self::OperatorEnum m_Operator;

  /// Downsampling radius in X direction.
  int m_RadiusX;

  /// Downsampling radius in Y direction.
  int m_RadiusY;

  /// Downsampling radius in Z direction.
  int m_RadiusZ;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationRegionFilter_h_included_
