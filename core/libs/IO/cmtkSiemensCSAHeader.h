/*
//
//  Copyright 2004-2012 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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
//  $Revision: 4527 $
//
//  $LastChangedDate: 2012-10-02 10:16:03 -0700 (Tue, 02 Oct 2012) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkSiemensCSAHeader_h_included_
#define __cmtkSiemensCSAHeader_h_included_

#include <cmtkconfig.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>

namespace
cmtk
{

/// Class for handling Siemens CSA headers in DICOM files.
class SiemensCSAHeader : public std::map< std::string,std::vector<std::string> >
{
public:
  /// This class.
  typedef SiemensCSAHeader Self;

  /// Constructor from binary blob.
  SiemensCSAHeader( const char* csaData, const size_t csaLength );
};

} // namespace cmtk

/// Write header contents to stream.
std::ostream& operator<<( std::ostream& stream, const cmtk::SiemensCSAHeader& csaHeader );
  

#endif // #ifndef __cmtkSiemensCSAHeader_h_included_
