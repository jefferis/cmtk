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

#ifndef __cmtkImageListEntry_h_included_
#define __cmtkImageListEntry_h_included_

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Class for entries of the image list.
   *@name ImageListEntry
   */
class ImageListEntry : 
  public Object
{
protected:
  /// Name of the current image.
  char *Name;

  /** Table (z) position of the image with respect to the select scan 
   * direction.
   */
  double Position;

  /** Default constructor.
   * Creates an empty entry and should therefore not be used. Consequently,
   * it's declared protected.
   */
  ImageListEntry() { Name = NULL; Position = 0; }
    
public:
  /** Create constructor.
   * Given a null-terminated string and a table position, this constructor
   * creates a matching image list entry.
   *@param name The image name. This pointer is freed upon destruction of
   * the list entry object. It must therefore have been allocated using
   * alloc() or strdup() and may not be used after creating this entry 
   * object.
   */
  ImageListEntry( char *name, const double position = 0 ) 
  { 
    Name = name; 
    Position = position; 
  }

  /** Copy operator.
   * String "Name" is duplicated by a call to strdup().
   */
  ImageListEntry& operator=( const ImageListEntry& other ) 
  {
    if ( other.Name ) 
      Name = strdup( other.Name );
    else
      Name = NULL;
    Position = other.Position;
    return *this;
  }

protected:
  /** Destructor.
   * If this entry's image name is non-NULL, the string is freed.
   */
  ~ImageListEntry() { if ( Name ) free(Name); }

  friend class Study;
};

//@}

} // namespace cmtk
  
#endif // #ifndef __cmtkImageListEntry_h_included_
