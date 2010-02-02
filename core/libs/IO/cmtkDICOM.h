/*
//
//  Copyright 2004-2010 SRI International
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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkDICOM_h_included_
#define __cmtkDICOM_h_included_

#include <cmtkconfig.h>

#include <cmtkImageIO.h>

#include <cmtkImageInfo.h>

#include <cmtkStudy.h>
#include <cmtkScalarImage.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Reader/writer class for DICOM images.
 */
class DICOM : 
  /// Inherit image IO interface.
  public ImageIO 
{
public:
  /// Return flags indicating which information is natively.
  virtual byte GetFormatCapabilities() const 
  { 
    return IMAGEFORMAT_DIMS | IMAGEFORMAT_STRUCTURE | IMAGEFORMAT_CALIBRATION;
  }

  /** Read ScalarImage from DICOM file.
   *@see ImageIO#Read
   */
  virtual ScalarImage* Read( const char *path, const Study* study = NULL, const int index = 0 ) const;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDICOM_h_included_
