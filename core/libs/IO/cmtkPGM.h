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

#ifndef __cmtkPGM_h_included_
#define __cmtkPGM_h_included_

#include <cmtkconfig.h>

#include <cmtkImageIO.h>

#include <cmtkImageInfo.h>

#include <cmtkScalarImage.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Reader/writer class for PGM (portable graymap) files.
 */
class PGM : 
  /// Inherit interface for image IO.
  public ImageIO 
{
public:
  /// Return flags indicating which information is natively.
  virtual byte GetFormatCapabilities() const 
  {
    return IMAGEFORMAT_DIMS | IMAGEFORMAT_STRUCTURE;
  }
  
  /** Read image from PGM file directly.
   */
  static ScalarImage* Read ( const char* );

  /** Read image from RAW file.
   *@see ImageIO#Read
   */
  virtual ScalarImage* Read ( const char*, const Study*, const int = 0 ) const;

  /** Write image to PGM file.
   *@see ImageIO#Write
   */
  virtual void Write ( const char*, const ImageInfo&, const int = 0 );

  /** Write scalar image to PGM file.
   */
  static void Write( const char* filename, const ScalarImage *image, const Types::DataItem greyFrom, const Types::DataItem greyTo );

  /** Write scalar image to 16bit PGM file.
   */
  static void Write( const char* filename, const ScalarImage *image );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkPGM_h_included_
