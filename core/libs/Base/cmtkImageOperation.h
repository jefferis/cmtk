/*
//
//  Copyright 2009 SRI International
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
//  $Revision: 1150 $
//
//  $LastChangedDate: 2010-01-18 12:41:38 -0800 (Mon, 18 Jan 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkImageOperation_h_included_
#define __cmtkImageOperation_h_included_

#include <cmtkconfig.h>

#include <cmtkSmartPtr.h>
#include <cmtkUniformVolume.h>

#include <list>

namespace
cmtk
{

/** Image operation base class.
 * Classes derived from this base class are used to implement an ordered sequence of operations
 * primarily for the "convertx" command line tool.
 *
 *\warning This class is not thread-safe in the sense that the costructed operation sequence is
 * stored in a static member field of this class, m_ImageOperationList.
 */
class ImageOperation
{
public:
  /// This class.
  typedef ImageOperation Self;

  /// Smart pointer.
  typedef cmtk::SmartPointer<Self> SmartPtr;
  
  /// Apply this operation to an image in place.
  virtual cmtk::UniformVolume::SmartPtr Apply( cmtk::UniformVolume::SmartPtr& volume ) 
  {
    return volume;
  }

  /// Apply all operations in list.
  static cmtk::UniformVolume::SmartPtr ApplyAll( cmtk::UniformVolume::SmartPtr& volume ) 
  {
    for ( std::list<Self::SmartPtr>::iterator opIt = Self::m_ImageOperationList.begin(); opIt != Self::m_ImageOperationList.end(); ++opIt )
      {
      volume = (*opIt)->Apply( volume );
      }
    return volume;
  }
  
protected:
  /// List of image operations.
  static std::list<Self::SmartPtr> m_ImageOperationList;
};

} // namespace cmtk

#endif // #ifndef __cmtkImageOperation_h_included_
