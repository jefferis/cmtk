/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkStudy_h_included_
#define __cmtkStudy_h_included_

#include <cmtkconfig.h>

#include <cmtkMacros.h>

#include <cmtkVolume.h>
#include <cmtkUniformVolume.h>
#include <cmtkAffineXform.h>
#include <cmtkWarpXform.h>
#include <cmtkLandmarkList.h>
#include <cmtkFileFormat.h>

#include <cmtkSegmentationLabel.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Constants to identify color table.
 */
enum
{
/// Grayscale color table.
  PALETTE_GRAY = 0,
/// Red color table.
  PALETTE_RED = 1,
/// Green color table.
  PALETTE_GREEN = 2,
/// Blue color table.
  PALETTE_BLUE = 3,
/// Rainbow color table.
  PALETTE_RAINBOW = 4,
/// Color table for labels.
  PALETTE_LABELS = 5
};

/** Class for parameters of a general imaging study.
 */
class Study
{
  /// Path of this study in the file system.
  cmtkGetSetMacroString(FileSystemPath);

  /// Short, memorable name assigned to this study.
  cmtkGetSetMacroString(Name);

  /// Textual description of study file type.
  cmtkGetSetMacroString(Description);

  /// Textual description of study file type.
  cmtkGetSetMacroString(Modality);

  /// Volume data associated with this study.
  cmtkGetSetMacro(UniformVolume::SmartPtr,Volume);

  /// Landmark list.
  cmtkGetSetMacro(LandmarkList::SmartPtr,LandmarkList)

  /// Voxel dimensions of the volume image.
  DataGrid::IndexType m_Dims;

  /// Voxel dimensions of the volume image.
  cmtkGetSetMacro3Array(Types::Coordinate,Calibration);

  /// Flag for custom calibration.
  cmtkGetSetMacro(bool,CustomCalibration);

  /// Minimum value.
  cmtkGetSetMacro(Types::DataItem,MinimumValue);

  /// Maximum value.
  cmtkGetSetMacro(Types::DataItem,MaximumValue);

  /// Pixel padding value.
  cmtkGetSetMacro(bool,Padding);

  /// Pixel padding value.
  cmtkGetSetMacro(Types::DataItem,PaddingValue);

  /// Flag for user-defined colormap.
  cmtkGetSetMacro(bool,HaveUserColorMap);

  /// Index of colormap.
  cmtkGetSetMacro(char,StandardColormap);

  /// Is colormap reversed?
  cmtkGetSetMacro(bool,ReverseColormap);

  /// Value corresponding to "black".
  cmtkGetSetMacro(Types::DataItem,Black);

  /// Value corresponding to "white".
  cmtkGetSetMacro(Types::DataItem,White);

  /// Gamma value.
  cmtkGetSetMacro(double,Gamma);

  /// Index of currently displayed image.
  cmtkGetSetMacro(unsigned int,DisplayedImageIndex);

  /// Displayed image zoom.
  cmtkGetSetMacro(unsigned int,ZoomFactor);

  /// Slice normal coordinate axis.
  cmtkGetSetMacro(int,SliceNormal);

public:
  /// Smart pointer to Study.
  typedef SmartPointer<Study> SmartPtr;

  /// Default constructor.
  Study();

  /// Constructor: Construct study from image file.
  Study( const char* fileSystemPath, const char* name = NULL );

  /// Destructor.
  virtual ~Study();

  /// Update from volume data, possibly after the data has been changed.
  virtual void UpdateFromVolume();

  /** Read volume data.
   *@param reRead If this is false, then the volume is only read if it has not
   * been read before. Otherwise, it is re-read in any case.
   *@return True if reading was successful; the "Volume" field has a pointer to
   * the resulting image volume.
   */
  virtual bool ReadVolume( const bool reRead = false, const char* orientation = NULL );

  /** Set study name; create name if no name given.
   * This function sets the name of this study. If no name is given (name
   * parameter is NULL pointer), then a name is constructed from the file 
   * system path of this study.
   */
  const char* SetMakeName( const char* name = NULL, const int suffix = 0 );

  /// Static study reader function.
  static Study* Read( const char* path );

  /** Copy colormap information from another study object.
   */
  virtual void CopyColormap( const Study* other );

  /** Get colormap from label list.
   */
  void SetFromLabelMap( const SegmentationLabelMap& lblMap ) 
  {
    this->m_HaveUserColorMap = true;
    UserLabelMap = lblMap;
  }

  /// Get user-defined label map.
  const SegmentationLabelMap& GetUserLabelMap() const 
  {
    return UserLabelMap; 
  }
  
private:
  /// User-defined label and colormap.
  SegmentationLabelMap UserLabelMap;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkStudy_h_included_

