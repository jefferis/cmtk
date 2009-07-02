/*
 *
 *  Copyright (C) 2002-2005, OFFIS
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
 *  Module:  dcmdata
 *
 *  Author:  Marco Eichelberg
 *
 *  Purpose: DcmOutputFileStream and related classes,
 *    implements streamed output to files.
 *
 *  Last Update:      $Author: meichel $
 *  Update Date:      $Date: 2005/12/08 15:41:22 $
 *  CVS/RCS Revision: $Revision: 1.6 $
 *  Status:           $State: Exp $
 *
 *  CVS/RCS Log at end of file
 *
 */

#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dcostrmf.h"
#include "dcmtk/ofstd/ofconsol.h"
#include "dcmtk/dcmdata/dcerror.h"

#define INCLUDE_CSTDIO
#define INCLUDE_CERRNO
#include "dcmtk/ofstd/ofstdinc.h"


DcmFileConsumer::DcmFileConsumer(const char *filename)
: DcmConsumer()
, file_(NULL)
, status_(EC_Normal)
{
  file_ = fopen(filename, "wb");
  if (!file_)
  {
    const char *text = strerror(errno);
    if (text == NULL) text = "(unknown error code)";
    status_ = makeOFCondition(OFM_dcmdata, 19, OF_error, text);
  }
}

DcmFileConsumer::DcmFileConsumer(FILE *file)
: DcmConsumer()
, file_(file)
, status_(EC_Normal)
{
}

DcmFileConsumer::~DcmFileConsumer()
{
  if (file_) fclose(file_);
}

OFBool DcmFileConsumer::good() const
{
  return status_.good();
}

OFCondition DcmFileConsumer::status() const
{
  return status_;
}

OFBool DcmFileConsumer::isFlushed() const
{
  return OFTrue;
}

Uint32 DcmFileConsumer::avail() const
{
  return OFstatic_cast(Uint32, -1); // assume unlimited file size
}

Uint32 DcmFileConsumer::write(const void *buf, Uint32 buflen)
{
  Uint32 result = 0;
  if (status_.good() && file_ && buf && buflen)
  {
    result = OFstatic_cast(Uint32, fwrite(buf, 1, OFstatic_cast(size_t, buflen), file_));
  }
  return result;
}

void DcmFileConsumer::flush()
{
  // nothing to flush
}

/* ======================================================================= */

DcmOutputFileStream::DcmOutputFileStream(const char *filename)
: DcmOutputStream(&consumer_) // safe because DcmOutputStream only stores pointer
, consumer_(filename)
{
}

DcmOutputFileStream::DcmOutputFileStream(FILE *file)
: DcmOutputStream(&consumer_) // safe because DcmOutputStream only stores pointer
, consumer_(file)
{
}

DcmOutputFileStream::~DcmOutputFileStream()
{
  // last attempt to flush stream before file is closed
  flush();
#ifdef DEBUG
  if (! isFlushed())
  {
    ofConsole.lockCerr() << "Warning: closing unflushed DcmOutputFileStream, loss of data!" << endl;
    ofConsole.unlockCerr();
  }
#endif
}


/*
 * CVS/RCS Log:
 * $Log: dcostrmf.cc,v $
 * Revision 1.6  2005/12/08 15:41:22  meichel
 * Changed include path schema for all DCMTK header files
 *
 * Revision 1.5  2004/02/04 16:36:47  joergr
 * Adapted type casts to new-style typecast operators defined in ofcast.h.
 *
 * Revision 1.4  2003/11/07 13:49:09  meichel
 * Added constructor taking an open FILE* instead of a file name string
 *
 * Revision 1.3  2002/11/27 12:06:50  meichel
 * Adapted module dcmdata to use of new header file ofstdinc.h
 *
 * Revision 1.2  2002/09/19 08:32:28  joergr
 * Added explicit type casts to keep Sun CC 2.0.1 quiet.
 *
 * Revision 1.1  2002/08/27 16:55:53  meichel
 * Initial release of new DICOM I/O stream classes that add support for stream
 *   compression (deflated little endian explicit VR transfer syntax)
 *
 *
 */
