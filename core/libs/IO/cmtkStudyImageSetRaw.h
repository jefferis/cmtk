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

#ifndef __cmtkStudyImageSetRaw_h_included_
#define __cmtkStudyImageSetRaw_h_included_

#include <cmtkconfig.h>

#include <cmtkStudyImageSet.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Study object for image sets composed of raw images.
class StudyImageSetRaw :
  /// Inherit basic image set study functions.
  public StudyImageSet
{
private:
  /// Convenience typedef.
  typedef StudyImageSet Superclass;

public:
  /// Length of file header before actual image data.
  igsGetSetMacro(unsigned int,HeaderLength);

  /// Bytes per image pixel.
  igsGetSetMacro(unsigned int,BytesPerPixel);

  /// Bytes per image pixel.
  igsGetSetMacro(bool,BigEndian);

  /// Bytes per image pixel.
  igsGetSetMacro(bool,Signed);

  /// Default constructor.
  StudyImageSetRaw() :
    HeaderLength( 0 ),
    BytesPerPixel( 1 ),
    BigEndian( false ),
    Signed( false )
  {}
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkStudyImageSetRaw_h_included_
