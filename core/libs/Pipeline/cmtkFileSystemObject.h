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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkFileSystemObject_h_included_
#define __cmtkFileSystemObject_h_included_

#include <cmtkPipelineObject.h>

#include <stdlib.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Class for objects with an associated filesystem representation.
 * This class provides information about whether the respective object has been
 * written to the filesystem, whether that representation is up-to-date, etc.
 */
class FileSystemObject : 
  public PipelineObject 
{
public:
  /// Return virtual class name.
  virtual const char* GetClassName() const
  { 
    return "cmtkFileSystemObject"; 
  }

  /** Return saved status of this object.
   */
  int GetIsSaved() const 
  {
    return SavedTime > this->GetModifiedTime();
  }

  /** Set saved status of this object.
   */
  void SetIsSaved() 
  {
    SavedTime = this->GetCurrentTime();
  }

  /** Get name of the filesystem represenation of this object.
   */
  const char* GetFileSystemName() const 
  {
    return FileSystemName;
  }

  /** Set name of the filesystem representation.
   */
  void SetFileSystemName( const char *fileSystemName );

  /** Is this object in its initial state?
   */
  int GetIsInitial() const 
  {
    return InitTime >= this->GetModifiedTime();
  }

  /** Notify initialization of this object.
   */
  void SetInitial() 
  {
    InitTime = this->GetCurrentTime();
  }
  
protected:
  /// Default constructor.
  FileSystemObject() 
  {
    SavedTime = -1;
    FileSystemName = NULL;
  }

  /// Virtual destructor.
  virtual ~FileSystemObject() 
  {
    if ( FileSystemName ) 
      free( FileSystemName );
  }

private:
  /** Time of the last save-to-disk operation.
   * Comparing this value to the last modification time allows to determine
   * whether this object has an up-to-date filesystem copy.
   */
  long SavedTime;

  /** Time when this object was last initialized.
   */
  long InitTime;
  
  /** Name of the filesystem represantation of this object.
   */
  char *FileSystemName;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFileSystemObject_h_included_
