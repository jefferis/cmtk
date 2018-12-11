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

#ifndef __cmtkImageOperationMapValues_h_included_
#define __cmtkImageOperationMapValues_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#include <map>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: apply mapping function to replace image pixel values.
class ImageOperationMapValues
  : public ImageOperation
{
public:
  /// This class.
  typedef ImageOperationMapValues Self;

  /// Superclass.
  typedef ImageOperation Superclass;

  /// Constructor.
  ImageOperationMapValues( const char* mapping /*!< Mapping function defined as 'VAL0[,VAL1,...][:NEWVAL]' to map values VAL0, VAL1, etc. to new value NEWVAL. If NEWVAL is not given, values are set to padding. */,
			   const bool exclusive = false /*!< Exclusive mapping flag: if set, all pixels not explicitly mapped will be set to padding; if not set, such pixels will keep their original values. */ );
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr  Apply( cmtk::UniformVolume::SmartPtr& volume );
  
  /// Create new operation to replace padded pixels.
  static void New( const char* mapping /*!< Mapping function defined as 'VAL0[,VAL1,...][:NEWVAL]' to map values VAL0, VAL1, etc. to new value NEWVAL. If NEWVAL is not given, values are set to padding.*/ )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new Self( mapping ) ) );
  }
  
  /// Create new operation to replace padded pixels and set all unmapped pixels to padding.
  static void NewExclusive( const char* mapping /*!< Mapping function defined as 'VAL0[,VAL1,...][:NEWVAL]' to map values VAL0, VAL1, etc. to new value NEWVAL. If NEWVAL is not given, values are set to padding.*/ )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new Self( mapping, true /*exclusive*/ ) ) );
  }
  
private:
  /// Mapping.
  std::map<Types::DataItem,Types::DataItem> m_Mapping;

  /** Exclusive mapping flag.
    * If set, all pixels not explicitly mapped will be set to padding; if not set, such pixels will keep their original values.
    */
  bool m_Exclusive;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationMapValues_h_included_
