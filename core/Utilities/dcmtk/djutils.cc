/*
 *
 *  Copyright (C) 1997-2005, OFFIS
 *
 *  This software and supporting documentation were developed by
 *
 *    Kuratorium OFFIS e.V.
 *    Healthcare Information and Communication Systems
 *    Escherweg 2
 *    D-26121 Oldenburg, Germany
 *
 *  THIS SOFTWARE IS MADE AVAILABLE,  AS IS,  AND OFFIS MAKES NO  WARRANTY
 *  REGARDING  THE  SOFTWARE,  ITS  PERFORMANCE,  ITS  MERCHANTABILITY  OR
 *  FITNESS FOR ANY PARTICULAR USE, FREEDOM FROM ANY COMPUTER DISEASES  OR
 *  ITS CONFORMITY TO ANY SPECIFICATION. THE ENTIRE RISK AS TO QUALITY AND
 *  PERFORMANCE OF THE SOFTWARE IS WITH THE USER.
 *
 *  Module:  dcmjpeg
 *
 *  Author:  Norbert Olges, Marco Eichelberg
 *
 *  Purpose: (STATUS: OK)
 *
 *  Last Update:      $Author: meichel $
 *  Update Date:      $Date: 2005/12/08 15:43:51 $
 *  Source File:      $Source: /share/dicom/cvs-depot/dcmtk/dcmjpeg/libsrc/djutils.cc,v $
 *  CVS/RCS Revision: $Revision: 1.4 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */

#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmjpeg/djutils.h"
#include "dcmtk/dcmdata/dcdeftag.h"  /* for tag constants */
#include "dcmtk/dcmdata/dcitem.h"    /* for class DcmItem */

#define INCLUDE_CCTYPE
#include "dcmtk/ofstd/ofstdinc.h"

const OFConditionConst EJC_Suspension(                           OFM_dcmjpeg,  1, OF_error, "IJG codec suspension return"  );
const OFConditionConst EJC_IJG8_FrameBufferTooSmall(             OFM_dcmjpeg,  2, OF_error, "Buffer for decompressed image (8 bits/sample) too small"  );
const OFConditionConst EJC_IJG12_FrameBufferTooSmall(            OFM_dcmjpeg,  3, OF_error, "Buffer for decompressed image (12 bits/sample) too small"  );
const OFConditionConst EJC_IJG16_FrameBufferTooSmall(            OFM_dcmjpeg,  4, OF_error, "Buffer for decompressed image (16 bits/sample) too small"  );
const OFConditionConst EJC_UnsupportedPhotometricInterpretation( OFM_dcmjpeg,  5, OF_error, "Codec does not support this PhotometricInterpretation"  );
const OFConditionConst EJC_UnsupportedColorConversion(           OFM_dcmjpeg,  6, OF_error, "Codec does not support this kind of color conversion"  );

const OFCondition      EJ_Suspension(                           EJC_Suspension);
const OFCondition      EJ_IJG8_FrameBufferTooSmall(             EJC_IJG8_FrameBufferTooSmall);
const OFCondition      EJ_IJG12_FrameBufferTooSmall(            EJC_IJG12_FrameBufferTooSmall);
const OFCondition      EJ_IJG16_FrameBufferTooSmall(            EJC_IJG16_FrameBufferTooSmall);
const OFCondition      EJ_UnsupportedPhotometricInterpretation( EJC_UnsupportedPhotometricInterpretation);
const OFCondition      EJ_UnsupportedColorConversion(           EJC_UnsupportedColorConversion);
EP_Interpretation DcmJpegHelper::getPhotometricInterpretation(DcmItem *item)
{
  if (item)
  {
    OFString photometric;
    if ((item->findAndGetOFString(DCM_PhotometricInterpretation, photometric)).good() && (photometric.length() > 0))
    {
      const char *c = photometric.c_str(); // guaranteed to be zero-terminated
      char cp[17]; // legal CS cannot be larger than 16 characters plus 0 byte
      int i=0; // current character
      while (*c && (i<16))
      {
        if (isalpha(*c)) cp[i++] = toupper(*c);
        else if (isdigit(*c)) cp[i++] = *c;
        c++;
      }
      cp[i] = 0; // write terminating zero
      photometric = cp; // copy back into OFString

      // now browse PhotometricInterpretationNames
      i = 0;
      while (PhotometricInterpretationNames[i].Name)
      {
        if (photometric == PhotometricInterpretationNames[i].Name) return PhotometricInterpretationNames[i].Type;
        i++;
      }
    }
  }
  return EPI_Unknown;
}


/*
 * CVS/RCS Log
 * $Log: djutils.cc,v $
 * Revision 1.4  2005/12/08 15:43:51  meichel
 * Changed include path schema for all DCMTK header files
 *
 * Revision 1.3  2005/11/30 14:13:14  onken
 * Added OFCondition constant for "unsupported color space conversions"
 *
 * Revision 1.2  2002/11/27 15:40:01  meichel
 * Adapted module dcmjpeg to use of new header file ofstdinc.h
 *
 * Revision 1.1  2001/11/13 15:58:35  meichel
 * Initial release of module dcmjpeg
 *
 *
 */
