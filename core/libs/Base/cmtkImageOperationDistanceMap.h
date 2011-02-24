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

#ifndef __cmtkImageOperationDistanceMap_h_included_
#define __cmtkImageOperationDistanceMap_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>
#include <Base/cmtkUniformDistanceMap.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Compute distance map.
class ImageOperationDistanceMap
/// Inherit generic image operation.
  : public ImageOperation
{
public:
  /// Distance map type.
  typedef UniformDistanceMap<double> DistanceMapType;

  /// Constructor.
  ImageOperationDistanceMap( const bool signedDistance ) : m_SignedDistance( signedDistance ) {}

  /// Apply this operation to an image in place.
  virtual UniformVolume::SmartPtr Apply( UniformVolume::SmartPtr& volume );

  /// Create new signed distance map operation.
  static void NewSigned()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationDistanceMap( true /*signedDistance*/) ) );
  }

  /// Create new unsigned distance map operation.
  static void NewUnsigned()
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationDistanceMap( false /*signedDistance*/) ) );
  }
  
private:
  /// Flag for signed (inside/outside) vs. unsigned (outside only) distance map.
  bool m_SignedDistance;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationDistanceMap_h_included_
