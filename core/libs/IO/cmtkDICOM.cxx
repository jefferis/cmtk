/*
//
//  Copyright 2014-2015 Google Inc.
//
//  Copyright 2004-2014 SRI International
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

#include "cmtkDICOM.h"

#include <Base/cmtkSurfaceNormal.h>
#include <Base/cmtkTypedArray.h>

#include <System/cmtkConsole.h>
#include <System/cmtkStrUtility.h>

#include <IO/cmtkSiemensCSAHeader.h>

#ifdef CMTK_USE_DCMTK_JPEG
#include <dcmtk/dcmjpeg/djdecode.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#include <stdio.h>
#include <string.h>
#include <ctime>

#ifndef DCM_ACR_NEMA_ImagePosition
#define DCM_ACR_NEMA_ImagePosition DcmTagKey(0x0020, 0x0030)
#endif

#ifndef DCM_ACR_NEMA_ImageOrientation
#define DCM_ACR_NEMA_ImageOrientation DcmTagKey(0x0020, 0x0035)
#endif

#ifndef DCM_ACR_NEMA_Location
#define DCM_ACR_NEMA_Location DcmTagKey(0x0020, 0x0050)
#endif

#ifndef DCM_ACR_NEMA_2C_VariablePixelData
#define DCM_ACR_NEMA_2C_VariablePixelData DcmTagKey(0x7f00, 0x0010)
#endif

namespace cmtk {

/** \addtogroup IO */
//@{

void DICOM::InitFromFile(const std::string &path) {
  this->m_Path = path;

#ifdef CMTK_USE_DCMTK_JPEG
  // register global decompression codecs
  static bool decodersRegistered = false;
  if (!decodersRegistered) {
    DJDecoderRegistration::registerCodecs(EDC_photometricInterpretation,
                                          EUC_default, EPC_default, 1);
    decodersRegistered = true;
  }
#endif

  std::auto_ptr<DcmFileFormat> fileformat(new DcmFileFormat);
  if (!fileformat.get()) {
    throw Exception("Could not create DICOM file format object.");
  }

  OFCondition status = fileformat->loadFile(path.c_str());

  if (!status.good()) {
    throw Exception("Cannot read DICOM file..");
    //    StdErr << "Error: cannot read DICOM file " << path << " (" <<
    //    status.text() << ")\n";
  }

  this->m_Dataset = fileformat->getAndRemoveDataset();

  if (!this->m_Dataset) {
    throw Exception("File format has NULL dataset.");
  }

  this->m_Document = std::auto_ptr<DiDocument>(
      new DiDocument(this->m_Dataset, this->m_Dataset->getOriginalXfer(),
                     CIF_AcrNemaCompatibility));
  if (!this->m_Document.get() || !this->m_Document->good()) {
    throw Exception("Could not create document representation.");
  }
}

const FixedVector<3, int> DICOM::GetDims() const {
  FixedVector<3, int> dims(0);

  Uint16 tempUint16 = 1;
  if (this->Document().getValue(DCM_Columns, tempUint16)) {
    dims[0] = static_cast<int>(tempUint16);
  }

  if (this->Document().getValue(DCM_Rows, tempUint16)) {
    dims[1] = static_cast<int>(tempUint16);
  }

  // detect and treat multi-frame files
  if (!this->Document().getValue(DCM_NumberOfFrames, tempUint16)) {
    // unlike Rows/Columns, NumberofFrames defaults to 1
    tempUint16 = 1;
  }
  dims[2] = tempUint16;

  return dims;
}

const FixedVector<3, double> DICOM::GetPixelSize() const {
  FixedVector<3, double> pixelSize(0.0);

  // get calibration from image
  const bool hasPixelSpacing =
      (this->Document().getValue(DCM_PixelSpacing, pixelSize[0], 0) > 0);
  if (hasPixelSpacing) {
    if (this->Document().getValue(DCM_PixelSpacing, pixelSize[1], 1) < 2) {
      throw Exception(
          "DICOM file does not have two elements in pixel size tag");
    }
  } else
    throw Exception("DICOM file does not specify pixel size");

  // get slice spacing from multi-slice images.
  if (!this->Document().getValue(DCM_SpacingBetweenSlices, pixelSize[2])) {
    pixelSize[2] = 0;
  }

  return pixelSize;
}

const FixedVector<3, double> DICOM::GetImageOrigin() const {
  FixedVector<3, double> imageOrigin(0.0);

  const char *image_position_s = NULL;
  if (!this->Document().getValue(DCM_ImagePositionPatient, image_position_s)) {
    // ImagePositionPatient tag not present, try ImagePosition instead
#ifdef DCM_ImagePosition
    if (!this->Document().getValue(DCM_ImagePosition, image_position_s))
      image_position_s = NULL;
#else
    if (!this->Document().getValue(DCM_ACR_NEMA_ImagePosition,
                                   image_position_s))
      image_position_s = NULL;
#endif
  }

  if (image_position_s) {
    double xyz[3];
    if (3 == sscanf(image_position_s, "%20lf%*c%20lf%*c%20lf", xyz, xyz + 1,
                    xyz + 2)) {
      imageOrigin = FixedVector<3, double>::FromPointer(xyz);
    }
  }

  return imageOrigin;
}

const FixedArray<2, FixedVector<3, double>> DICOM::GetImageOrientation() const {
  FixedArray<2, FixedVector<3, double>> orientation;

  orientation[0] = FixedVector<3, double>(0.0);
  orientation[1] = FixedVector<3, double>(0.0);

  orientation[0][0] = 1;
  orientation[1][1] = 1;

  const char *image_orientation_s = NULL;
#ifdef DCM_ImageOrientation
  if (!this->Document().getValue(DCM_ImageOrientation, image_orientation_s))
#else
  if (!this->Document().getValue(DCM_ACR_NEMA_ImageOrientation,
                                 image_orientation_s))
#endif
  {
    // ImageOrientation tag not present, try ImageOrientationPatient instead
    if (!this->Document().getValue(DCM_ImageOrientationPatient,
                                   image_orientation_s))
      image_orientation_s = NULL;
  }

  if (image_orientation_s) {
    double dx[3], dy[3];
    if (6 == sscanf(image_orientation_s,
                    "%20lf%*c%20lf%*c%20lf%*c%20lf%*c%20lf%*c%20lf", dx, dx + 1,
                    dx + 2, dy, dy + 1, dy + 2)) {
      orientation[0] = (FixedVector<3, double>::FromPointer(dx));
      orientation[1] = (FixedVector<3, double>::FromPointer(dy));
    }
  }

  return orientation;
}

TypedArray::SmartPtr DICOM::GetPixelDataArray(const size_t pixelDataLength) {
  DcmElement *delem = NULL;

  unsigned short bitsAllocated = 0;
  if ((delem = this->m_Document->search(DCM_BitsAllocated))) {
    delem->getUint16(bitsAllocated);
  } else {
    // No "BitsAllocated" tag; use "BitsStored" instead.
    if ((delem = this->m_Document->search(DCM_BitsStored))) {
      delem->getUint16(bitsAllocated);
    }
  }

  bool pixelDataSigned = false;
  Uint16 pixelRepresentation = 0;
  if (this->m_Document->getValue(DCM_PixelRepresentation, pixelRepresentation) >
      0)
    pixelDataSigned = (pixelRepresentation == 1);

  double rescaleIntercept, rescaleSlope;
  const bool haveRescaleIntercept =
      (0 != this->m_Document->getValue(DCM_RescaleIntercept, rescaleIntercept));
  if (!haveRescaleIntercept) rescaleIntercept = 0;

  const bool haveRescaleSlope =
      (0 != this->m_Document->getValue(DCM_RescaleSlope, rescaleSlope));
  if (!haveRescaleSlope) rescaleSlope = 1;

  pixelDataSigned = pixelDataSigned || (rescaleIntercept < 0);

  Uint16 paddingValue = 0;
  const bool paddingFlag =
      (this->m_Dataset->findAndGetUint16(DCM_PixelPaddingValue, paddingValue))
          .good();

  TypedArray::SmartPtr pixelDataArray;

#ifdef DCM_VariablePixelData
  delem = this->m_Document->search(DCM_VariablePixelData);
#else
  delem = this->m_Document->search(DCM_ACR_NEMA_2C_VariablePixelData);
#endif
  if (!delem) delem = this->m_Document->search(DCM_PixelData);

  if (delem) {
    if ((delem->getTag().getEVR() == EVR_OW) || (bitsAllocated > 8)) {
      Uint16 *pdata = NULL;
      delem->getUint16Array(pdata);
      if (pixelDataSigned) {
        const short paddingShort = static_cast<short>(paddingValue);
        pixelDataArray = TypedArray::Create(
            TYPE_SHORT, pdata, pixelDataLength, paddingFlag, &paddingShort,
            Memory::ArrayCXX::DeleteWrapper<short>);
      } else {
        const unsigned short paddingUShort =
            static_cast<unsigned short>(paddingValue);
        pixelDataArray = TypedArray::Create(
            TYPE_USHORT, pdata, pixelDataLength, paddingFlag, &paddingUShort,
            Memory::ArrayCXX::DeleteWrapper<unsigned short>);
      }
    } else {
      Uint8 *pdata = NULL;
      delem->getUint8Array(pdata);
      if (pixelDataSigned) {
        const char paddingChar = static_cast<char>(paddingValue);
        pixelDataArray = TypedArray::Create(
            TYPE_CHAR, pdata, pixelDataLength, paddingFlag, &paddingChar,
            Memory::ArrayCXX::DeleteWrapper<char>);
      } else {
        const char paddingByte = static_cast<byte>(paddingValue);
        pixelDataArray = TypedArray::Create(
            TYPE_BYTE, pdata, pixelDataLength, paddingFlag, &paddingByte,
            Memory::ArrayCXX::DeleteWrapper<byte>);
      }
    }

    delem->detachValueField();
  }

  if (!pixelDataArray) {
    throw Exception("Could not read pixel data from DICOM file");
  }

  if (haveRescaleIntercept || haveRescaleSlope) {
    double intpart = 0;
    if (fabs(modf(rescaleSlope, &intpart) / rescaleSlope) > 1e-5) {
      pixelDataArray = pixelDataArray->Convert(TYPE_FLOAT);
    }

    pixelDataArray->Rescale(rescaleSlope, rescaleIntercept);
  }

  return pixelDataArray;
}

const FixedVector<3, double> DICOM::DemosaicAndGetNormal(
    const FixedArray<2, FixedVector<3, double>> &imageOrientation,
    const FixedVector<3, Types::Coordinate> &deltas, FixedVector<3, int> &dims,
    TypedArray::SmartPtr &pixelDataArray, FixedVector<3, double> &imageOrigin) {
  // without further information, we "guess" the image normal vector
  FixedVector<3, double> sliceNormal =
      SurfaceNormal(imageOrientation[0], imageOrientation[1]).Get();

  // detect and treat Siemens multi-slice mosaics
  const char *tmpStr = NULL;
  if (this->Document().getValue(DCM_Manufacturer, tmpStr)) {
    if (!strncmp(tmpStr, "SIEMENS", 7)) {
      Uint16 tempUint16 = 0;
      const DcmTagKey nSlicesTag(0x0019, 0x100a);
      if (this->Document().getValue(nSlicesTag, tempUint16)) {
        dims[2] = tempUint16;
      }

      // check for mosaic
      if (dims[2] || (this->Document().getValue(DCM_ImageType, tmpStr) &&
                      strstr(tmpStr, "MOSAIC"))) {
        int unmosaicImageRows;
        int unmosaicImageCols;

        const DcmTagKey mosaicTag(0x0051, 0x100b);
        if (this->Document().getValue(mosaicTag, tmpStr)) {
          if (2 != sscanf(tmpStr, "%6dp*%6ds", &unmosaicImageRows,
                          &unmosaicImageCols)) {
            if (2 != sscanf(tmpStr, "%6d*%6ds", &unmosaicImageRows,
                            &unmosaicImageCols)) {
              StdErr
                  << "ERROR: unable to parse mosaic size from (0x0051,0x100b): "
                  << tmpStr << "\n";
            }
          }
        }

        // For the following, see here:
        // http://nipy.sourceforge.net/nibabel/dicom/siemens_csa.html#csa-header
        this->ParseSiemensCSA(DcmTagKey(0x0029, 0x1020), unmosaicImageCols,
                              unmosaicImageRows, dims[2], sliceNormal,
                              imageOrigin);  // series information
        this->ParseSiemensCSA(DcmTagKey(0x0029, 0x1010), unmosaicImageCols,
                              unmosaicImageRows, dims[2], sliceNormal,
                              imageOrigin);  // image information

        // hopefully we have figured out the mosaic dimensions by now.
        if ((unmosaicImageCols > 0) && (unmosaicImageRows > 0)) {
          const int xMosaic = dims[0] / unmosaicImageCols;

          dims[0] = unmosaicImageCols;
          dims[1] = unmosaicImageRows;

          // de-mosaic the data array
          const unsigned long imageSizePixels = dims[0] * dims[1] * dims[2];
          TypedArray::SmartPtr newDataArray(
              TypedArray::Create(pixelDataArray->GetType(), imageSizePixels));

          const size_t pixelsPerSlice = unmosaicImageCols * unmosaicImageRows;
          size_t toOffset = 0;
          for (int slice = 0; slice < dims[2]; ++slice) {
            for (int j = 0; j < unmosaicImageRows; ++j, toOffset += dims[0]) {
              const size_t iPatch = slice % xMosaic;
              const size_t jPatch = slice / xMosaic;

              const size_t fromOffset = jPatch * xMosaic * pixelsPerSlice +
                                        j * xMosaic * unmosaicImageCols +
                                        iPatch * unmosaicImageCols;
              pixelDataArray->BlockCopy(*newDataArray, toOffset, fromOffset,
                                        unmosaicImageCols);
            }
          }

          pixelDataArray = newDataArray;

          // convert CSA-header center-of-image origin to corner-of-image
          // standard origin (Issue #6754)
          imageOrigin -=
              (0.5 * ((dims[0] - 1) * deltas[0] * imageOrientation[0] +
                      (dims[1] - 1) * deltas[1] * imageOrientation[1]));
        }
      }
    }
  }
  return sliceNormal;
}

void DICOM::ParseSiemensCSA(const DcmTagKey &tagKey, int &unmosaicImageCols,
                            int &unmosaicImageRows, int &slices,
                            FixedVector<3, double> &sliceNormal,
                            FixedVector<3, double> &imageOrigin) {
  const Uint8 *csaHeaderInfo = NULL;
  unsigned long csaHeaderLength = 0;
  if (this->Dataset()
          .findAndGetUint8Array(tagKey, csaHeaderInfo, &csaHeaderLength)
          .status() == OF_ok) {
    SiemensCSAHeader csaHeader((const char *)csaHeaderInfo, csaHeaderLength);

    SiemensCSAHeader::const_iterator it =
        csaHeader.find("AcquisitionMatrixText");
    if ((it != csaHeader.end()) && !it->second.empty()) {
      if (2 != sscanf(it->second[0].c_str(), "%6dp*%6ds", &unmosaicImageRows,
                      &unmosaicImageCols)) {
        if (2 != sscanf(it->second[0].c_str(), "%6d*%6ds", &unmosaicImageRows,
                        &unmosaicImageCols)) {
          StdErr << "ERROR: unable to parse mosaic size from CSA field "
                    "AcquisitionMatrixText: "
                 << it->second[0] << " in file " << this->m_Path << "\n";
        }
      }
    }

    it = csaHeader.find("NumberOfImagesInMosaic");
    if ((it != csaHeader.end()) && !it->second.empty())
      slices = atof(it->second[0].c_str());

    it = csaHeader.find("SliceNormalVector");
    if ((it != csaHeader.end()) && (it->second.size() >= 3)) {
      for (size_t i = 0; i < 3; ++i)
        sliceNormal[i] = atof(it->second[i].c_str());
    }

    // get true slice0 location from CSA header
    it = csaHeader.find("MrPhoenixProtocol");
    if ((it != csaHeader.end()) && !it->second.empty()) {
      // ID strings for three axes in CSA header
      const std::string sliceOrientationString[] = {"dSag", "dCor", "dTra"};

      for (int i = 0; i < 3; ++i) {
        const size_t sliceTagPos = it->second[0].find(
            std::string("sSliceArray.asSlice[0].sPosition.") +
            sliceOrientationString[i]);
        if (sliceTagPos != std::string::npos) {
          const size_t equalPos = it->second[0].find('=', sliceTagPos);
          if (equalPos != std::string::npos) {
            imageOrigin[i] = atof(it->second[0].substr(equalPos + 1).c_str());
          } else {
            StdErr << "ERROR: unable to get image origin component from: "
                   << it->second[0] << " in file " << this->m_Path
                   << "\nAssuming zero.\n";
            imageOrigin[i] = 0;
          }
        } else {
          StdErr << "ERROR: unable to get image origin tag for component "
                 << sliceOrientationString[i] << " from CSA header in file "
                 << this->m_Path << "\nAssuming zero.\n";
          imageOrigin[i] = 0;
        }
      }
    }
  }
}

ScalarImage *DICOM::Read(const char *path) {
  ScalarImage *image = NULL;

  Self dicom(path);

  FixedVector<3, int> dims = dicom.GetDims();
  FixedVector<3, double> pixelSize = dicom.GetPixelSize();
  ScalarImage::SpaceVectorType imageOrigin = dicom.GetImageOrigin();

  image = new ScalarImage(dims[0], dims[1], dims[2]);
  image->SetPixelSize(pixelSize[0], pixelSize[1]);
  image->SetFrameToFrameSpacing(pixelSize[2]);

  TypedArray::SmartPtr pixelDataArray =
      dicom.GetPixelDataArray(dims[0] * dims[1] * dims[2]);
  image->SetPixelData(pixelDataArray);

  // now some more manual readings...

  // get original table position from image.
  double sliceLocation = 0;
  if (!dicom.Document().getValue(DCM_SliceLocation, sliceLocation)) {
#ifdef DCM_Location
    dicom.Document().getValue(DCM_Location, sliceLocation);
#else
    dicom.Document().getValue(DCM_ACR_NEMA_Location, sliceLocation);
#endif
  }
  image->SetImageSlicePosition(sliceLocation);
  image->SetImageOrigin(imageOrigin);

  // get original image direction from file.
  FixedArray<2, FixedVector<3, double>> imageOrientation =
      dicom.GetImageOrientation();
  image->SetImageDirectionX(imageOrientation[0]);
  image->SetImageDirectionY(imageOrientation[1]);

  return image;
}

//@}

}  // namespace cmtk
