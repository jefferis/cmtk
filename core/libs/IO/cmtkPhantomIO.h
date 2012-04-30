/*
//
//  Copyright 2012 SRI International
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

#ifndef __cmtkPhantomIO_h_included_
#define __cmtkPhantomIO_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkDetectedPhantomMagphanEMR051.h>

#include <mxml.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Read and write imaging phantom descriptions to and from XML files.
class PhantomIO
{
public:
  /// This class.
  typedef PhantomIO Self;
  
  /// Write Magphan EMR051 description.
  static void Write( const DetectedPhantomMagphanEMR051& phantom, const std::string& fpath );

private:
  /// Whitespace callback function for MiniXML.
  static const char* WhitespaceWriteMiniXML( mxml_node_t* node, int where);
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkPhantomIO_h_included_
