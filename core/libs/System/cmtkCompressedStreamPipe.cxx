/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include "cmtkCompressedStream.h"

#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

namespace cmtk {

/** \addtogroup System */
//@{

CompressedStream::Pipe::Pipe(const std::string &filename, const char *command) {
  char cmd[PATH_MAX];

#ifndef _MSC_VER
  if (static_cast<size_t>(snprintf(cmd, sizeof(cmd), command,
                                   filename.c_str())) >= sizeof(cmd)) {
    StdErr
        << "WARNING: length of path exceeds system PATH_MAX in "
           "CompressedStream::OpenDecompressionPipe and will be truncated.\n";
  }
  errno = 0;

  this->m_File = popen(cmd, CMTK_FILE_MODE);
  if (!this->m_File) {
    fprintf(stderr, "ERROR: popen(\"%s\") returned NULL (errno=%d).\n", cmd,
            errno);
    perror("System message");
    throw 0;
  }
#else
  this->m_TempName[0] = 0;
  if (snprintf(cmd, sizeof(cmd), command, filename, tmpnam(this->m_TempName)) >=
      sizeof(cmd)) {
    StdErr
        << "WARNING: length of path exceeds system PATH_MAX in "
           "CompressedStream::OpenDecompressionPipe and will be truncated.\n";
  }

  _flushall();
  int sysReturn = system(cmd);

  if (sysReturn) {
    fprintf(stderr, "Command %s returned %d\n", cmd, sysReturn);
    fprintf(stderr, "Errno = %d\n", errno);
  }

  this->m_File = fopen(this->m_TempName, CMTK_FILE_MODE);

  if (!this->m_File) {
    throw 0;
  }
#endif

  this->m_BytesRead = 0;
}

CompressedStream::Pipe::~Pipe() { this->Close(); }

void CompressedStream::Pipe::Close() {
#ifndef _MSC_VER
  if (this->m_File) {
    pclose(this->m_File);
    this->m_File = NULL;
  }
#else
  if (this->m_TempName[0]) {
    remove(this->m_TempName);
    this->m_TempName[0] = 0;
  }
#endif  // # ifndef _MSC_VER
}

void CompressedStream::Pipe::Rewind() {
  StdErr << "CompressedStream::Pipe::Rewind() is not implemented\n";
  throw ExitException(1);
}

size_t CompressedStream::Pipe::Read(void *data, size_t size, size_t count) {
  const size_t result = fread(data, size, count, this->m_File);
  this->m_BytesRead += result;
  return result / size;
}

bool CompressedStream::Pipe::Get(char &c) {
  const int data = fgetc(this->m_File);
  if (data != EOF) {
    c = (char)data;
    ++this->m_BytesRead;
    return true;
  }

  return false;
}

int CompressedStream::Pipe::Tell() const { return this->m_BytesRead; }

bool CompressedStream::Pipe::Feof() const { return (feof(this->m_File) != 0); }

}  // namespace cmtk
