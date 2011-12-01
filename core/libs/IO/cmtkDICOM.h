/*
//
//  Copyright 2004-2011 SRI International
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

#include <Base/cmtkScalarImage.h>

#include <dcmtk/dcmdata/dcdeftag.h>
#include <dcmtk/dcmimgle/didocu.h>
#include <dcmtk/dcmimgle/diutils.h>

#include <memory>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Reader/writer class for DICOM images.
 */
class DICOM
{
public:
  /// This class.
  typedef DICOM Self;

  /** Constructor.
   */
  DICOM( const char* path );

  /// Get image dimensions (number of pixels per axis).
  const FixedVector<3,int> GetDims() const;

  /// Get image pixel size.
  const FixedVector<3,double> GetPixelSize() const;

  /** Get pixel data array.
   * The pixel data type is determined automatically based on bits allocated and signed vs. unsigned representation.
   * 
   * If the RescaleSlope or RescaleIntercept tags are present, intensity rescaling is applied.
   *
   * If a padding value is defined in the DICOM file, this value is also set as padding in the output array.
   *
   *\warning As a side effect, this function releases the pixel data array pointer from the DICOM object, i.e., 
   *  this function can only be called ONCE for each object.
   */
  TypedArray::SmartPtr GetPixelDataArray( const size_t pixelDataLength );

  /// Get const DICOM dataset.
  const DcmDataset& Dataset() const
  {
    return *(this->m_Dataset);
  }

  /// Get DICOM dataset.
  DcmDataset& Dataset()
  {
    return *(this->m_Dataset);
  }

  /// Get const DICOM document.
  const DiDocument& Document() const
  {
    return *(this->m_Document);
  }

  /// Get DICOM document.
  DiDocument& Document()
  {
    return *(this->m_Document);
  }

  /** Read ScalarImage from DICOM file.
   */
  static ScalarImage* Read( const char *path );

private:
  /// Pointer to the DICOM dataset object
  DcmDataset* m_Dataset;

  /// Pointer to the DICOM document object
  std::auto_ptr<DiDocument> m_Document;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDICOM_h_included_
