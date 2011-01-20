/*
//
//  Copyright 2010 Torsten Rohlfing
//
//  Copyright 2011 SRI International
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

#ifndef __cmtkImageOperationReplace_h_included_
#define __cmtkImageOperationReplace_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

namespace
cmtk
{

/// Image operation: replace image pixel values.
class ImageOperationReplace
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// This class.
  typedef ImageOperationReplace Self;

  /// Superclass.
  typedef ImageOperation Superclass;

  /// Operation mode.
  typedef enum
  {
    /// Replace padded pixels.
    REPLACE_PADDING,
    /// Replace Inf and NaN pixels.
    REPLACE_INF_NAN
  } Mode;

  /// Constructor.
  ImageOperationReplace( const Self::Mode mode /*!< Operation mode.*/, const Types::DataItem value /*!< Replacement data value.*/ )
    : m_Mode( mode ),
      m_ReplacementValue( value )
  {}

  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    TypedArray& volumeData = *(volume->GetData());
    switch ( this->m_Mode ) 
      {
      case Self::REPLACE_PADDING:
	volumeData.ReplacePaddingData( this->m_ReplacementValue );
	break;
      case Self::REPLACE_INF_NAN:
#pragma omp parallel for
	for ( size_t i = 0; i < volumeData.GetDataSize(); ++i )
	  {
	  cmtk::Types::DataItem value = 0;
	  if ( volumeData.Get( value, i ) )
	    {
	    if ( !finite( value ) )
	      {
	      volumeData.Set( this->m_ReplacementValue, i );
	      }
	    }
	  }
	break;
      }
    return volume;
  }
  
  /// Create new operation to replace padded pixels.
  static void NewReplacePadding( const double value /*!< Replacement value. */ )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new Self( Self::REPLACE_PADDING, value ) ) );
  }

  /// Create new operation to replace pixels with Inf or NaN values.
  static void NewReplaceInfNaN( const double value /*!< Replacement value. */ )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new Self( Self::REPLACE_INF_NAN, value ) ) );
  }
  
private:
  /// Operation mode.
  Self::Mode m_Mode;

  /// Padding value.
  Types::DataItem m_ReplacementValue;
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationReplace_h_included_
