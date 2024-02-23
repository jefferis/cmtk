/*
//
//  Copyright 2004-2014, 2022 SRI International
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

#include "cmtkImageFileDICOM.h"

#include <System/cmtkDebugOutput.h>

#include <IO/cmtkFileFormat.h>

#include <dcmtk/dcmimgle/diutils.h>

#include <sstream>
#include <string>
#include <numeric>
#include <algorithm>

namespace cmtk
{

  void
  ImageFileDICOM::Print() const
  {
    cmtk::DebugOutput(1) << "  File Name =            [" << this->m_FileDir << "/" << this->m_FileName << "]\n";
    cmtk::DebugOutput(1) << "  SeriesID =             [" << this->GetTagValue(DCM_SeriesInstanceUID) << "]\n";
    cmtk::DebugOutput(1) << "  StudyID =              [" << this->GetTagValue(DCM_StudyInstanceUID) << "]\n";
    cmtk::DebugOutput(1) << "  ImagePositionPatient = [" << this->GetTagValue(DCM_ImagePositionPatient) << "]\n";
    cmtk::DebugOutput(1) << "  AcquisitionNumber =    [" << this->m_AcquisitionNumber << "]\n";
    cmtk::DebugOutput(1) << "  Modality =             [" << this->GetTagValue(DCM_Modality) << "]\n";

    if (this->GetTagValue(DCM_Modality) == "MR")
    {
      cmtk::DebugOutput(1) << "  EchoTime =          [" << this->GetTagValue(DCM_EchoTime) << "]\n";
      cmtk::DebugOutput(1) << "  RepetitionTime =      [" << this->GetTagValue(DCM_RepetitionTime) << "]\n";
    }
  }

  bool
  ImageFileDICOM::Match(const Self &other, const Types::Coordinate numericalTolerance, const bool disableCheckOrientation, const bool ignoreAcquisitionNumber) const
  {
    // do not stack multislice images
    if (this->m_IsMultislice || other.m_IsMultislice)
      return false;

    if (!disableCheckOrientation)
    {
      double orientThis[6];
      if (6 != sscanf(this->GetTagValue(DCM_ImageOrientationPatient).c_str(), "%30lf%*c%30lf%*c%30lf%*c%30lf%*c%30lf%*c%30lf", orientThis, orientThis + 1, orientThis + 2, orientThis + 3, orientThis + 4, orientThis + 5))
      {
        StdErr << "ERROR: unable to parse image orientation from '" << this->GetTagValue(DCM_ImageOrientationPatient).c_str() << "'\n";
        return false;
      }

      double orientOther[6];
      if (6 != sscanf(other.GetTagValue(DCM_ImageOrientationPatient).c_str(), "%30lf%*c%30lf%*c%30lf%*c%30lf%*c%30lf%*c%30lf", orientOther, orientOther + 1, orientOther + 2, orientOther + 3, orientOther + 4, orientOther + 5))
      {
        StdErr << "ERROR: unable to parse image orientation from '" << other.GetTagValue(DCM_ImageOrientationPatient).c_str() << "'\n";
        return false;
      }

      for (int i = 0; i < 6; ++i)
      {
        if (fabs(orientThis[i] - orientOther[i]) > numericalTolerance)
          return false;
      }
    }

    return (this->m_FileDir == other.m_FileDir) && // do not stack across directories
           (this->GetTagValue(DCM_FrameOfReferenceUID) == other.GetTagValue(DCM_FrameOfReferenceUID)) &&
           (this->GetTagValue(DCM_SeriesInstanceUID) == other.GetTagValue(DCM_SeriesInstanceUID)) &&
           (this->GetTagValue(DCM_StudyInstanceUID) == other.GetTagValue(DCM_StudyInstanceUID)) &&
           (this->GetTagValue(DCM_EchoTime) == other.GetTagValue(DCM_EchoTime)) &&
           (this->GetTagValue(DCM_RepetitionTime) == other.GetTagValue(DCM_RepetitionTime)) &&
           (this->GetTagValue(DCM_InversionTime) == other.GetTagValue(DCM_InversionTime)) &&
           (this->m_BValue == other.m_BValue) &&
           (this->m_BVector == other.m_BVector) &&
           ((this->m_AcquisitionNumber == other.m_AcquisitionNumber) || ignoreAcquisitionNumber) &&
           (this->m_RawDataType == other.m_RawDataType);
  }

  ImageFileDICOM::ImageFileDICOM(const std::string &filepath)
      : m_IsMultislice(false),
        m_IsDWI(false),
        m_DwellTime(0.0),
        m_PhaseEncodeDirectionSign(""),
        m_BValue(0),
        m_BVector(0.0),
        m_HasBVector(false),
        m_RawDataType("unknown")
  {
    if (cmtk::FileFormat::Identify(filepath, false /*decompress*/) != cmtk::FILEFORMAT_DICOM) // need to disable "decompress" in Identify() because DCMTK cannot currently read using on-the-fly decompression.
      throw(0);

    this->m_FileName = filepath;
    this->m_FileDir = "";

    const size_t lastSlash = this->m_FileName.rfind(CMTK_PATH_SEPARATOR);
    if (lastSlash != std::string::npos)
    {
      this->m_FileDir = this->m_FileName.substr(0, lastSlash);

      // trim trailing slashes, both forward and back
      const size_t lastNotSlash = this->m_FileDir.find_last_not_of("/\\");
      if (lastNotSlash != std::string::npos)
      {
        this->m_FileDir.erase(lastNotSlash + 1);
      }
      else
      {
        this->m_FileDir.clear();
      }

      this->m_FileName = this->m_FileName.substr(lastSlash + 1);

      const size_t suffix = this->m_FileName.rfind('.');
      if (suffix != std::string::npos)
      {
        const std::string suffixStr = this->m_FileName.substr(suffix + 1);
        if ((suffixStr == ".Z") || (suffixStr == ".gz"))
        {
          this->m_FileName = this->m_FileName.erase(suffix);
        }
      }
    }

    std::auto_ptr<DcmFileFormat> fileformat(new DcmFileFormat);
    OFCondition status = fileformat->loadFile(filepath.c_str());

    if (!status.good())
    {
      cmtk::StdErr << "Error: cannot read DICOM file " << filepath << " (" << status.text() << ")\n";
      throw(0);
    }

    this->m_Dataset = std::auto_ptr<DcmDataset>(fileformat->getAndRemoveDataset());
    if (!this->m_Dataset.get())
    {
      throw(1);
    }

    this->m_Document = std::auto_ptr<DiDocument>(new DiDocument(this->m_Dataset.get(), this->m_Dataset->getOriginalXfer(), CIF_AcrNemaCompatibility));
    if (!this->m_Document.get() || !this->m_Document->good())
    {
      throw(2);
    }

    // read the string tags that we need for later
    const DcmTagKey defaultStringTags[] = {
        DCM_PixelSpacing,
        DCM_Manufacturer, DCM_ManufacturerModelName, DCM_DeviceSerialNumber, DCM_StationName, DCM_SoftwareVersions, DCM_AcquisitionTime,
        DCM_PatientsName,
        DCM_Modality, DCM_EchoTime, DCM_RepetitionTime, DCM_InversionTime, DCM_ImagingFrequency, DCM_SequenceName, DCM_InPlanePhaseEncodingDirection,
        DCM_FlipAngle, DCM_MagneticFieldStrength, DCM_PhotometricInterpretation, DCM_SliceThickness,
        DCM_StudyInstanceUID, DCM_StudyID, DCM_StudyDate,
        DCM_FrameOfReferenceUID, DCM_SeriesInstanceUID, DCM_SeriesDescription,
        DCM_ImagePositionPatient, DCM_ImageOrientationPatient,
        DCM_RescaleIntercept, DCM_RescaleSlope,
        // GE-specific tags
        DCM_GE_PulseSequenceName, DCM_GE_PulseSequenceDate, DCM_GE_InternalPulseSequenceName, DCM_GE_AssetRFactors,
        // Siemens-specific tags

        // Philips-specific tags

        // Zero-tag ends the array
        DcmTagKey(0, 0)};

    for (size_t tagIdx = 0; (defaultStringTags[tagIdx].getGroup() != 0) && (defaultStringTags[tagIdx].getElement() != 0); ++tagIdx)
    {
      const char *tmpStr = NULL;
      if (this->m_Document->getValue(defaultStringTags[tagIdx], tmpStr))
        this->m_TagToStringMap[defaultStringTags[tagIdx]] = tmpStr;
    }

    this->DoFOV();

    // read acquisition matrix
    this->DoAcquisitionMatrix();

    // check for multi-slice DICOMs
    Uint16 nFrames = 0;
    if (this->m_Document->getValue(DCM_NumberOfFrames, nFrames))
    {
      this->m_IsMultislice = (nFrames > 1);
    }

    // read integer tags that we also need
    if (!this->m_Document->getValue(DCM_InstanceNumber, this->m_InstanceNumber))
      this->m_InstanceNumber = 0;

    if (!this->m_Document->getValue(DCM_AcquisitionNumber, this->m_AcquisitionNumber))
      this->m_AcquisitionNumber = 0;

    // check for which vendor and deal with specifics elsewhere
    std::string vendor( this->m_TagToStringMap[DCM_Manufacturer] );
    std::transform(vendor.begin(), vendor.end(), vendor.begin(), ::toupper);
    if ( vendor.compare( 0, 7, "SIEMENS" ) == 0 )
    {
      this->DoVendorTagsSiemens();
    }
    else if ( vendor.compare( 0, 7, "PHILIPS" ) == 0 )
    {
      this->DoVendorTagsPhilips();
    }  
    else if ( vendor.compare( 0, 2, "GE" ) == 0 )
    {
      this->DoVendorTagsGE();
    }
  }

  void
  ImageFileDICOM::DoFOV()
  {
    const DcmTagKey intTagKeys[] = {
      DCM_Columns, DCM_Rows,

      DcmTagKey(0, 0)
    };

    for (size_t tagIdx = 0; (intTagKeys[tagIdx].getGroup() != 0) && (intTagKeys[tagIdx].getElement() != 0); ++tagIdx)
    {
      const DcmTagKey& tagKey = intTagKeys[tagIdx];
      const Uint16 *pTmpInt = NULL;
      DcmObject *pObj = NULL;
      if (this->m_Document->getValue(tagKey, pTmpInt, pObj) > 0) {
        std::ostringstream s;
        s << *pTmpInt;
        this->m_TagToStringMap[tagKey] = std::string(s.str());
      }
    }
  }

  template<typename InputIt>
  std::string join(InputIt first, InputIt last, const std::string& separator = ",")
  {
    std::ostringstream result;
    if (first != last) {
      result << *first;
      while (++first != last) {
        result << separator << *first;
      }
    }
    return result.str();
  }

  void
  ImageFileDICOM::DoAcquisitionMatrix()
  {
    Uint16 tmpDbl = 0;
    DcmObject *pObj = NULL;
    std::vector<std::string> vAcquisitionMatrix;
    for (size_t dim = 0, nDims = this->m_Document->getVM(DCM_AcquisitionMatrix); dim < nDims; ++dim) 
    {
      if (this->m_Document->getValue(DCM_AcquisitionMatrix, tmpDbl, dim, pObj) > 0) 
      {
        std::ostringstream s;
        s << tmpDbl;
        vAcquisitionMatrix.push_back(std::string(s.str()));
      }
    }

    this->m_TagToStringMap[DCM_AcquisitionMatrix] = join(vAcquisitionMatrix.begin(), vAcquisitionMatrix.end(), "\\");
  }

  bool 
  ImageFileDICOM::DoVendorTagsSiemensDwellTime(const SiemensCSAHeader* csaImageHeader)
  {
    // Get the dwell time from the dedicated "RealDwellTime" private  
    // SIEMENS DICOM attribute if it exists, or from the CSA header if it  
    // exists and contains it. The units of the SIEMENS dwell times are
    // nano-seconds.
    bool dwell_time_found = false;
    double dwell_time_value;
    unsigned long dwell_time_length = 0;
    const Uint8* dwell_time_buffer = NULL;
    OFString dwell_time_str;

    if ( this->m_Dataset->findAndGetUint8Array( DCM_Siemens_RealDwellTime, dwell_time_buffer, &dwell_time_length ).good() ) 
    {
      if( dwell_time_found = ((dwell_time_buffer != NULL) && (dwell_time_length > 0)) )
      {
        dwell_time_str.assign( OFstatic_cast(const char*, dwell_time_buffer), dwell_time_length);
        dwell_time_value = atof( dwell_time_str.c_str() );
      }
    } 

    if ( !dwell_time_found && this->m_Document->getValue(DCM_Siemens_RealDwellTime, dwell_time_str) > 0 )
    {
      if ( dwell_time_found = !dwell_time_str.empty() ) 
        dwell_time_value = atof( dwell_time_str.c_str() );
    }

    if ( !dwell_time_found && csaImageHeader != NULL )
    {
      SiemensCSAHeader::const_iterator csa_it = csaImageHeader->find( "RealDwellTime" );
      if ( dwell_time_found = ((csa_it != csaImageHeader->end()) && !csa_it->second.empty()) )
        dwell_time_value = atof( csa_it->second[0].c_str() );
    }
    // I do not know why the reciprocal of the reported dwell time is being
    // stored in the dwell time member. Does this not give us Bandwidth
    // instead of dwell time?
    if ( dwell_time_found ) 
    {
      if ( dwell_time_value != 0.0 ) 
        this->m_DwellTime = 1.0 / dwell_time_value;
      else
        this->m_DwellTime = 0.0;
    }
    return dwell_time_found;
  }

  bool 
  ImageFileDICOM::DoVendorTagsSiemensPhaseEncodeDirection(const SiemensCSAHeader* csaImageHeader)
  {
    // Get the polarity of the phase encoding from the dedicated 
    // "PhaseEncodingDirectionPositive" private SIEMENS DICOM attribute if 
    // it exists, or from the CSA header if it exists and contains it
    bool polarity_found = false;
    char polarity_value = '\0';
    unsigned long polarity_length = 0;
    const Uint8* polarity_buffer = NULL;
    OFString polarity_str;

    if ( this->m_Dataset->findAndGetUint8Array( DCM_Siemens_PhaseEncodingDirectionPositive, polarity_buffer, &polarity_length ).good() ) 
    {
      if( polarity_found = (polarity_buffer != NULL && polarity_length > 0) )
        polarity_value = *( OFstatic_cast(const char*, polarity_buffer) );
    }

    if ( !polarity_found && this->m_Document->getValue(DCM_Siemens_PhaseEncodingDirectionPositive, polarity_str) )
    {
      if ( polarity_found = !polarity_str.empty() )
        polarity_value = polarity_str[0];
    }

    if ( !polarity_found && csaImageHeader != NULL )
    {
      SiemensCSAHeader::const_iterator csa_it = csaImageHeader->find( "PhaseEncodingDirectionPositive");
      if ( polarity_found = (csa_it != csaImageHeader->end() && !csa_it->second.empty()) )
        polarity_value = csa_it->second[0][0];

    }
    // use of the polarity_found flag here is to preserve the previous
    // functionality where the m_PhaseEncodeDirectionSign member was
    // only set when polarity_value was found. 
    if ( polarity_found )
    {
      switch ( polarity_value )
      {
      case '0':
        this->m_PhaseEncodeDirectionSign = "NEG";
        break;
      case '1':
        this->m_PhaseEncodeDirectionSign = "POS";
        break;
      default:
        this->m_PhaseEncodeDirectionSign = "UNKNOWN";	    
        break;
      }
    }
    return polarity_found;
  }

  bool 
  ImageFileDICOM::DoVendorTagsSiemensDiffusionBValue(const SiemensCSAHeader* csaImageHeader)
  {
    unsigned long bvalue_length = 0;
    const Uint8* bvalue_buffer = NULL;
    OFString bvalue_str;

    if ( this->m_Dataset->findAndGetUint8Array( DCM_Siemens_DiffusionBValue, bvalue_buffer, &bvalue_length ).good() ) 
    {
      if( bvalue_buffer != NULL && bvalue_length > 0 )
      {
        bvalue_str.assign( OFstatic_cast(const char*, bvalue_buffer), bvalue_length);
        this->m_BValue = atof( bvalue_str.c_str() );
        return true;
      }
    }

    if ( this->m_Document->getValue(DCM_Siemens_DiffusionBValue, bvalue_str) > 0 )
    {
      if( !bvalue_str.empty() )
      {
        this->m_BValue = atof( bvalue_str.c_str() );
        return true;
      }
    }

    if ( csaImageHeader != NULL )
    {
      SiemensCSAHeader::const_iterator csa_it = csaImageHeader->find( "B_value" );
      if ( (csa_it != csaImageHeader->end()) && !csa_it->second.empty() )
      {
        this->m_BValue = atof( csa_it->second[0].c_str() );
        return true;
      }
    }
    return false;
  }

  bool 
  ImageFileDICOM::DoVendorTagsSiemensDiffusionBVector(const SiemensCSAHeader* csaImageHeader)
  {
    unsigned long bvector_length = 0;
    const Uint8* bvector_buffer = NULL;

    if ( this->m_Dataset->findAndGetUint8Array( DCM_Siemens_DiffusionGradientOrientation, bvector_buffer, &bvector_length ).good() ) 
    {
      if ( bvector_buffer != NULL && (bvector_length * sizeof(Uint8)) >= (this->m_BVector.Size() * sizeof(double)) )
      {
        const double* input_itr = OFstatic_cast(const double*, bvector_buffer);
        std::copy(input_itr, input_itr + this->m_BVector.Size(), this->m_BVector.begin());
        return true;
      }  
    }

    if ( this->m_Document->getValue(DCM_Siemens_DiffusionGradientOrientation, this->m_BVector[0], 0) >= this->m_BVector.Size() )
    {
      for ( size_t idx = 1; idx < this->m_BVector.Size(); ++idx )
        this->m_Document->getValue(DCM_Siemens_DiffusionGradientOrientation, this->m_BVector[idx], idx);
      return true;
    }

    if ( csaImageHeader != NULL )
    {
      SiemensCSAHeader::const_iterator csa_it = csaImageHeader->find( "DiffusionGradientDirection" );	
      if ( (csa_it != csaImageHeader->end()) && (csa_it->second.size() >= this->m_BVector.Size()) )
      {
        for ( size_t idx = 0; idx < this->m_BVector.Size(); ++idx )
          this->m_BVector[idx] = atof( csa_it->second[idx].c_str() );
        return true;
      }
    }
    return false;
  }

  void
  ImageFileDICOM::DoVendorTagsSiemens()
  {
    Uint16 nFrames = 0;
    double tmpDbl = 0;
    const char *tmpStr = NULL;
    SiemensCSAHeader::const_iterator csa_it;

    this->m_IsMultislice = (0 != this->m_Document->getValue(DcmTagKey(0x0019, 0x100a), nFrames));            // Number of Slices tag
    this->m_IsMultislice |= (this->m_Document->getValue(DCM_ImageType, tmpStr) && strstr(tmpStr, "MOSAIC")); // mosaics are always multi-slice

    // Retrieve slice times if this tag exists (it's an array of double-precision floats)
    for (size_t slice = 0, nslices = this->m_Document->getVM(DCM_Siemens_MosaicRefAcqTimes); slice < nslices; ++slice)
    {
      if (this->m_Document->getValue(DCM_Siemens_MosaicRefAcqTimes, tmpDbl, slice) > 0)
        this->m_SliceTimes.push_back(tmpDbl);
    }

    if (this->GetTagValue(DCM_Modality) == "MR")
    {
      // try to extract raw data type from "ImageType"
      if (this->m_Document->getValue(DCM_ImageType, tmpStr))
      {
        if (strstr(tmpStr, "\\P\\"))
          this->m_RawDataType = "phase";
        else if (strstr(tmpStr, "\\M\\"))
          this->m_RawDataType = "magnitude";
        else if (strstr(tmpStr, "\\R\\"))
          this->m_RawDataType = "real";
      }

      // If the CSA Header exists, instantiate SiemensCSAHeader to use in
      // the retrieval of subsequent imagining quantities
      const Uint8* csaImageHeaderInfo = NULL;
      unsigned long csaImageHeaderLength = 0;
      std::auto_ptr<SiemensCSAHeader> csaImageHeader;
      if ( this->m_Dataset->findAndGetUint8Array( DcmTagKey(0x0029,0x1010), csaImageHeaderInfo, &csaImageHeaderLength ).status() == OF_ok ) // the "Image" CSA header, not the "Series" header.
      {
        csaImageHeader.reset( new SiemensCSAHeader( OFstatic_cast(const char*, csaImageHeaderInfo), csaImageHeaderLength ) );
      }

      DoVendorTagsSiemensDwellTime(csaImageHeader.get());
      DoVendorTagsSiemensPhaseEncodeDirection(csaImageHeader.get());
 
      bool bvalue_found = DoVendorTagsSiemensDiffusionBValue(csaImageHeader.get());
      // 20161028 djk: make m_IsDWI is always true if the B_value field is defined
      this->m_IsDWI |= bvalue_found;

      this->m_HasBVector = DoVendorTagsSiemensDiffusionBVector(csaImageHeader.get());
      this->m_HasBVector |= bvalue_found;
      this->m_IsDWI |= this->m_HasBVector;

      bool found_directionality = false;
      OFString directionality;
      unsigned long directionality_length = 0;
      const Uint8* directionality_buffer = NULL;
      if ( this->m_Dataset->findAndGetUint8Array( DCM_Siemens_DiffusionDirectionality, directionality_buffer, &directionality_length ).good() ) 
      {
        if( found_directionality = (directionality_buffer != NULL && directionality_length > 0) )
          directionality.assign( OFstatic_cast(const char *, directionality_buffer), directionality_length );
      }

      if ( !found_directionality )
        found_directionality = (this->m_Document->getValue(DCM_Siemens_DiffusionDirectionality, directionality) > 0);

      if ( !found_directionality  && csaImageHeader.get() != NULL ) 
      {
        csa_it = csaImageHeader->find( "DiffusionDirectionality" );
        if ( found_directionality = (csa_it != csaImageHeader->end() && !csa_it->second.empty()) )
          directionality = csa_it->second[0].c_str();
      }

      if ( found_directionality )
      {
        if ( 0 == directionality.compare( 0, 11, "DIRECTIONAL" ) ) 
        {
          this->m_IsDWI = true;
          this->m_HasBVector = true;
        }
        else if ( 0 == directionality.compare( 0, 9, "ISOTROPIC" ) ) 
        {
          this->m_IsDWI = true;
          this->m_HasBVector = false;
        }
      }
    }
  }

  void
  ImageFileDICOM::DoVendorTagsGE()
  {
    int tmpInt = 0;
    double tmpDbl = 0;

    if (this->GetTagValue(DCM_Modality) == "MR")
    {
      // raw data type
      Sint16 rawTypeIdx = 3;
      if (!this->m_Document->getValue(DCM_GE_RawDataType_ImageType, rawTypeIdx))
        rawTypeIdx = 0; // assume this is a magnitude image
      rawTypeIdx = std::min(3, std::max(0, (int)rawTypeIdx));

      const char *const RawDataTypeString[4] = {"magnitude", "phase", "real", "imaginary"};
      this->m_RawDataType = RawDataTypeString[rawTypeIdx];

      Sint16 effEchoSpacing = 0;
      if (this->m_Document->getValue(DCM_GE_EffectiveEchoSpacing, effEchoSpacing))
      {
        std::ostringstream ss;
        ss << effEchoSpacing;
        this->m_TagToStringMap[DCM_GE_EffectiveEchoSpacing] = ss.str();

        // Default (no acceleration) - dwell time is effective echo spacing (convert microseconds to seconds here)
        this->m_DwellTime = 1e-6 * effEchoSpacing;

        // Check for "asset" acceleration
        const std::string asset = this->GetTagValue(DCM_GE_AssetRFactors, "");
        if (asset != "")
        {
          float rFactor;
          if (1 == sscanf(asset.c_str(), "%10f\\%*c", &rFactor))
          {
            // With acceleration - dwell time is effective echo spacing multiplied by acceleration factor (e.g., 0.5 for 2-fold acceleration)
            this->m_DwellTime *= rFactor;
          }
        }
      }

      // dwi information
      this->m_IsDWI = false;
      const char *tmpStr = NULL;
      if (this->m_Document->getValue(DcmTagKey(0x0019, 0x10e0), tmpStr) > 0) // Number of Diffusion Directions
      {
        const int nDirections = atoi(tmpStr);
        if (nDirections > 0)
        {
          this->m_IsDWI = true;

          if (this->m_Document->getValue(DcmTagKey(0x0043, 0x1039), tmpStr) > 0) // bValue tag
          {
            if (1 == sscanf(tmpStr, "%10d\\%*c", &tmpInt))
            {
              this->m_BValue = static_cast<double>(tmpInt);

              this->m_HasBVector = true;
              for (int i = 0; i < 3; ++i)
              {
                if (this->m_Document->getValue(DcmTagKey(0x0019, 0x10bb + i), tmpStr) > 0) // bVector tags
                {
                  this->m_BVector[i] = atof(tmpStr);
                }
                else
                {
                  this->m_BVector[i] = 0;
                  this->m_HasBVector = false;
                }
              }
              this->m_BVector[2] *= -1; // for some reason z component apparently requires negation of sign on GE
            }
          }
        }
      }
    }
  }

  void
  ImageFileDICOM::DoVendorTagsPhilips()
  {
    double tmpDbl = 0;

    if (this->GetTagValue(DCM_Modality) == "MR")
    {
      if (this->m_Document->getValue(DcmTagKey(0x0018, 0x9087), tmpDbl) > 0) // bValue tag
      {
        this->m_IsDWI = true;
        this->m_BValue = tmpDbl;
      }

      this->m_HasBVector = true;
      if (this->m_BValue > 0)
      {
        for (size_t i = 0; this->m_IsDWI && (i < 3); ++i)
        {
          if (this->m_Document->getValue(DcmTagKey(0x0018, 0x9089), tmpDbl, i) > 0)
            this->m_BVector[i] = tmpDbl;
          else
            this->m_IsDWI = false;
        }

        const char *tmpStr = NULL;
        if (this->m_Document->getValue(DcmTagKey(0x2001, 0x1004), tmpStr) > 0) // DiffusionDirection tag - "O", oblique or "I", isotropic. Uses "I" with bValue!=0 to identify average diffusion image
        {
          if (tmpStr)
          {
            this->m_HasBVector = (tmpStr[0] != 'I');
          }
        }
      }
    }
  }

  bool
  ImageFileDICOM::MatchAnyPattern(const std::map<DcmTagKey, std::string> &patterns) const
  {
    // check for positive include list
    if (!patterns.empty())
    {
      for (std::map<DcmTagKey, std::string>::const_iterator it = patterns.begin(); it != patterns.end(); ++it)
      {
        // if tag not found, do not include
        const char *tmpStr = NULL;
        if (this->m_Document->getValue(it->first, tmpStr))
        {
          // if tag value matches, then return true
          if (strstr(tmpStr, it->second.c_str()))
            return true;
        }
      }
    }

    return false; // did not find a match or list was empty
  }

  bool
  ImageFileDICOM::MatchAllPatterns(const std::map<DcmTagKey, std::string> &patterns) const
  {
    // check for positive include list
    if (!patterns.empty())
    {
      for (std::map<DcmTagKey, std::string>::const_iterator it = patterns.begin(); it != patterns.end(); ++it)
      {
        // if tag not found, do not include
        const char *tmpStr = NULL;
        if (this->m_Document->getValue(it->first, tmpStr))
        {
          // if tag value matches, then return true
          if (!strstr(tmpStr, it->second.c_str()))
            return false;
        }
      }
    }

    return true; // all matched or list empty
  }

} // namespace CMTK
