/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkXformListIO_h_included__
#define __cmtkXformListIO_h_included__

#include <cmtkconfig.h>

#include <Base/cmtkXformList.h>

#include <vector>
#include <string>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Utility class to generate a list of concatenated transformation objects.
 */
class XformListIO 
{
public:
  /// This class.
  typedef XformListIO Self;

  /** Create transformation list from string list. 
   *\return An XformList object with concatenated Xform (and derived) objects, each of which
   * may be optionally inverted.
   */
  static XformList MakeFromStringList( const std::vector<std::string>& stringList /*!< List of transformation paths. If an entry is "--inverse" or "-i", then the next following transformation is marked to be applied inverse.*/ );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkXformListIO_h_included__
