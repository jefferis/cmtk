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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkStudyImageSet_h_included_
#define __cmtkStudyImageSet_h_included_

#include <cmtkconfig.h>

#include <cmtkStudy.h>

#include <list>
#include <string>

#include <cmtkFileFormat.h>

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
  igsGetSetMacro(bool,MultiFile);

  /// Directory that contains the image files.
  igsGetSetMacroString(ImageDirectory);

  /// Directory that contains the image files.
  igsGetSetMacro(FileFormatID,ImageFormat);

  /// Default constructor.
  StudyImageSet() : 
    Study(),
    MultiFile( false ),
    ImageDirectory( NULL ), 
    ImageFormat( FILEFORMAT_UNKNOWN )
  {}

  /** Read volume data.
   *@param reRead If this is false, then the volume is only read if it has not
   * been read before. Otherwise, it is re-read in any case.
   *@return True if reading was successful; the "Volume" field has a pointer to
   * the resulting image volume.
   *\note Parameter "orientation" is not used here; it is only present to
   *  avoid hiding the inherited virtual function.
   */
    virtual bool ReadVolume( const bool reRead = false, const char* = NULL );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkStudyImageSet_h_included_
