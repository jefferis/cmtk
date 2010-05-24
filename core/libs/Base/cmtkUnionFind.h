/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#ifndef __cmtkUnionFind_h_included_
#define __cmtkUnionFind_h_included_

#include <cmtkconfig.h>

#include <set>
#include <list>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Class template for (relatively) efficient uni-find algorithm.
template<class T>
class UnionFind
{
public:
  /// Internal set type.
  typedef std::set<T> SetType;

  /// Internal list type.
  typedef std::list<SetType> ListType;

  /// Opqaue result of "find" operation.
  typedef typename ListType::iterator FindResultType;

  /// Find operation.
  FindResultType Find( const T& key )
  {
    for ( FindResultType it = this->m_List.begin(); it != this->m_List.end(); ++it )
      {
      if ( it->find( key ) != it->end() )
	return it;
      }
    return this->End();
  }

  /// Find representative key.
  const T FindKey( const T& key )
  {
    return *(this->Find( key )->begin());
  }

  /// End-of-list iterator.
  FindResultType End()
  {
    return this->m_List.end();
  }

  /// Union operation.
  void Union( const FindResultType& s1, const FindResultType& s2 )
  {
    if ( s1 != s2 )
      {
      s1->insert( s2->begin(), s2->end() );
      this->m_List.erase( s2 );
      }
  }

  /// Insert a new key by itself.
  void Insert( const T& key )
  {
    SetType newSet;
    newSet.insert( key );
    this->m_List.push_back( newSet );
  }
  
private:
  /// The list of sets.
  ListType m_List;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkUnionFind_h_included_
