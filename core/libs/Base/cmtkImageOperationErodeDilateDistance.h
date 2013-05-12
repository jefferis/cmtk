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
//  $Revision: 2902 $
//
//  $LastChangedDate: 2011-02-24 12:12:46 -0800 (Thu, 24 Feb 2011) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkImageOperationErodeDilateDistance_h_included_
#define __cmtkImageOperationErodaDilateDistance_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageOperation.h>
#include <Base/cmtkUniformVolumeMorphologicalOperators.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Image operation: erode or dilate by distance (rather than pixels).
class ImageOperationErodeDilateDistance
/// Inherit from image operation base class.
  : public ImageOperation
{
public:
  /// Constructor:
  ImageOperationErodeDilateDistance( const double distance ) : m_Distance( distance ) {}
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr Apply( cmtk::UniformVolume::SmartPtr& volume )
  {
    if ( this->m_Distance < 0 )
      {
      cmtk::UniformVolumeMorphologicalOperators ops( volume );
      volume->SetData( ops.GetErodedByDistance( -this->m_Distance ) );
      }
    else
      {
      if ( this->m_Distance > 0 )
	{
	cmtk::UniformVolumeMorphologicalOperators ops( volume );
	volume->SetData( ops.GetDilatedByDistance( this->m_Distance ) );
	}
      }
    return volume;
  }

  /// Create new dilation operation.
  static void NewDilate( const double distance )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationErodeDilateDistance( distance ) ) );
  }

  /// Create new erosion operation.
  static void NewErode( const double distance )
  {
    ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationErodeDilateDistance( -distance ) ) );
  }
  
private:
  /// Distance of erosion (if negative) or dilation (if positive).
  double m_Distance;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageOperationErodeDilateDistance_h_included_
