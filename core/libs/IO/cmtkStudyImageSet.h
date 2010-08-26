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

#ifndef __cmtkStudyImageSet_h_included_
#define __cmtkStudyImageSet_h_included_

#include <cmtkconfig.h>

#include "IO/cmtkStudy.h"
#include "IO/cmtkFileFormat.h"

#include <list>
#include <string>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// An imaging study that is constructed from multiple 2-D images.
class StudyImageSet :
  /// Inherit basic study fields and functions.
  public Study,
  /// Inherit string list for list of file names.
  public std::list<std::string>
{
private:
  /// Convenience typedef.
  typedef Study Superclass;

public:
  /// Is this a single file or a multi-file study?
  cmtkGetSetMacro(bool,MultiFile);

  /// Directory that contains the image files.
  cmtkGetSetMacroString(ImageDirectory);

  /// Directory that contains the image files.
  cmtkGetSetMacro(FileFormatID,ImageFormat);

  /// Default constructor.
  StudyImageSet() : 
    Study(),
    m_MultiFile( false ),
    m_ImageDirectory( NULL ), 
    m_ImageFormat( FILEFORMAT_UNKNOWN )
  {}
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkStudyImageSet_h_included_
