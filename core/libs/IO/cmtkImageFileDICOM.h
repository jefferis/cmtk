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
//  $Revision: 4497 $
//
//  $LastChangedDate: 2012-08-24 13:46:21 -0700 (Fri, 24 Aug 2012) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkImageFileDICOM_h_included_
#define __cmtkImageFileDICOM_h_included_

#include <cmtkconfig.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkSmartConstPtr.h>

#include <Base/cmtkFixedVector.h>

#include <dcmtk/dcmdata/dctk.h>
#include <dcmtk/dcmimgle/didocu.h>

#include <string>
#include <map>
#include <memory>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Class handling a single DICOM image file and its meta data.
class ImageFileDICOM
{
public:
  /// This class.
  typedef ImageFileDICOM Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to constant object of this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// File name.
  char* fname;

  /// File system path (i.e., directory).
  char* fpath;

  /// Modality.
  std::string Modality;

  /// Device manufacturer.
  std::string Manufacturer;

  /// Device model.
  std::string ManufacturerModel;

  /// Device serial number.
  std::string m_DeviceSerialNumber;

  /// Flag for multislice images
  bool IsMultislice;

  /// Patient name.
  std::string PatientName;

  /// DICOM SeriesUID.
  std::string SeriesUID;

  /// DICOM FrameOfReferenceUID
  std::string m_FrameOfReferenceUID;

  /// DICOM SeriesDescription
  std::string SeriesDescription;

  /// DICOM StudyID.
  std::string StudyID;

  /// DICOM StudyDate.
  std::string StudyDate;
  
  /// DICOM StudyUID.
  std::string StudyUID;

  /// DICOM pixel value rescale intercept
  std::string m_RescaleIntercept;

  /// DICOM pixel value rescale slope
  std::string m_RescaleSlope;

  /// MR repetition time, TR
  std::string RepetitionTime;

  /// MR echo time, TE.
  std::string EchoTime;

  /// MR imaging frequency.
  std::string m_ImagingFrequency;

  /// 3D image position (first pixel) in patient coordinates.
  std::string ImagePositionPatient;

  /// 3D image orientation (pos. x and pos. y image axes) in patient coordinates.
  std::string ImageOrientationPatient;

  /// DICOM acquisition number.
  Sint32 AcquisitionNumber;

  /// DICOM image number (index in volume).
  Sint32 InstanceNumber;

  /// Raw data type (real, imaginary, phase, magnitude) currently supported on GE images only.
  std::string RawDataType;

  /// Flag for diffusion-weighted images.
  bool IsDWI;

  /// B value for DWI.
  Sint16 BValue;

  /// B vector for DWI.
  cmtk::FixedVector<3,double> BVector;

  /// Constructor.
  ImageFileDICOM( const char* filename );

  /// Destructor.
  ~ImageFileDICOM();

  /// Determine whether two images match, i.e., belong to the same volume.
  bool Match( const Self& other, const Types::Coordinate numericalTolerance = 0, /*!< Numerical comparison tolerance; values with absolute difference less than this threshold are considered equal. */
	      const bool disableCheckOrientation = false /*!< Flag for disabling the checking of image orientation vectors.*/,
	      const bool ignoreAcquisitionNumber = false /*!< When this flag is set, the AcquisitionNumber DICOM tag is ignore for matching images*/ ) const;

  /// Test if this image matches at least one from a list of DICOM tag patterns.
  bool MatchAnyPattern( const std::map<DcmTagKey,std::string>& patterns ) const;

  /// Test if this image matches all from a list of DICOM tag patterns (or list is empty).
  bool MatchAllPatterns( const std::map<DcmTagKey,std::string>& patterns ) const;

  /// Compare order based on file name (for lexicographic sorting).
  static bool lessFileName( const Self* lhs, const Self* rhs )
  {
    return strcmp( lhs->fname, rhs->fname ) < 0;
  }

  /// Compare order based on image instace (for sorting in acquisition order).
  static bool lessInstanceNumber( const Self* lhs, const Self* rhs )
  {
    return lhs->InstanceNumber < rhs->InstanceNumber;
  }

  /// Print informatiomn about this object.
  void Print() const;

  /// Release memory allocated for complete DICOM document.
  void ReleaseDocument()
  {
    this->m_Document = std::auto_ptr<DiDocument>( NULL );
  }

private:
  /// DICOM document object.
  std::auto_ptr<DiDocument> m_Document;

  /// Handle Siemens private tags.
  void DoVendorTagsSiemens();

  /// Handle GE private tags.
  void DoVendorTagsGE();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageFileDICOM_h_included_
