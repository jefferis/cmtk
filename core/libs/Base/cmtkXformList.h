/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkXformList_h_included_
#define __cmtkXformList_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkXformListEntry.h>

#include <System/cmtkSmartPtr.h>

#include <deque>

namespace
cmtk
{

/** \addtogroup Base */
//@{
/// A transformation list.
class XformList :
  /// Inherit STL list.
  public std::deque< XformListEntry::SmartConstPtr > 
{
private:
  /// Error threshold for inverse approximation.
  Types::Coordinate m_Epsilon;
  
public:
  /// This class.
  typedef XformList Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Constructor.
  XformList( const Types::Coordinate epsilon = 0.0 ) : m_Epsilon( epsilon ) {};
  
  /// Set epsilon.
  void SetEpsilon( const Types::Coordinate epsilon ) 
  {
    this->m_Epsilon = epsilon;
  }
  
  /// Add a transformation the the end of the list, i.e., to be applied after the current list of transformations
  void Add( const Xform::SmartConstPtr& xform, const bool inverse = false, const Types::Coordinate globalScale = 1.0 );
  
  /// Add a transformation the the end of the list, i.e., to be applied before the current list of transformations
  void AddToFront( const Xform::SmartConstPtr& xform, const bool inverse = false, const Types::Coordinate globalScale = 1.0 );
  
  /// Apply a sequence of (inverse) transformations.
  bool ApplyInPlace( Xform::SpaceVectorType& v ) const;
  
  /// Get the Jacobian determinant of a sequence of transformations.
  bool GetJacobian( const Xform::SpaceVectorType& v, Types::DataItem& jacobian, const bool correctGlobalScale = true ) const;

  /// Is this transformation list all affine?
  bool AllAffine() const;

  /// Make all-affine copy of this transformation list.
  Self MakeAllAffine() const;

  /** Get fixed image path, if available.
   * Not every transformation file format stores the fixed image path, in which case
   * an empty string is returned here.
   */
  std::string GetFixedImagePath() const;

  /** Get moving image path, if available.
   * Not every transformation file format stores the moving image path, in which case
   * an empty string is returned here.
   */
  std::string GetMovingImagePath() const;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkXformList_h_included_
