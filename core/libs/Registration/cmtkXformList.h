/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#ifdef HAVE_CMTKCONFIG_H
#  include <cmtkconfig.h>
#endif

#include <cmtkXformListEntry.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{
/// A transformation list.
class XformList :
  /// Inherit STL list.
  public std::list< SmartPointer<XformListEntry> > 
{
private:
  /// Error threshold for inverse approximation.
  Types::Coordinate m_Epsilon;
  
public:
  /// Constructor.
  XformList( const Types::Coordinate epsilon = 0.0 ) : m_Epsilon( epsilon ) {};
  
  /// Set epsilon.
  void SetEpsilon( const Types::Coordinate epsilon ) 
  {
    this->m_Epsilon = epsilon;
  }
  
  /// Add a transformation
  void Add( Xform::SmartPtr& xform, const bool inverse = false, const Types::Coordinate globalScale = 1.0 );
  
  /// Apply a sequence of (inverse) transformations.
  bool ApplyInPlace( Vector3D& v );
  
  /// Get the Jacobian determinant of a sequence of transformations.
  bool GetJacobian( const Vector3D& v, Types::DataItem& jacobian, const bool correctGlobalScale = true );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkXformList_h_included_
