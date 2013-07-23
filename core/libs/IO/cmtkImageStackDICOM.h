/*
//
//  Copyright 2004-2013 SRI International
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

#ifndef __cmtkImageStackDICOM_h_included_
#define __cmtkImageStackDICOM_h_included_

#include <cmtkconfig.h>

#include <IO/cmtkImageFileDICOM.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkSmartConstPtr.h>

#include <Base/cmtkUniformVolume.h>

#include <mxml.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Class handling a stack of DICOM image files.
class ImageStackDICOM : public std::vector<ImageFileDICOM::SmartConstPtr> 
{
public:
  /// This class.
  typedef ImageStackDICOM Self;
  
  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const object of this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Enum type to select DICOM information to be embedded into output images as "description".
  typedef enum
  {
    /// No embedding.
    EMBED_NONE = 0,
    /// Embed StudyID plus Date.
    EMBED_STUDYID_STUDYDATE = 1,
    /// Embed patient name.
    EMBED_PATIENTNAME = 2,
    /// Embed series description.
    EMBED_SERIESDESCR = 3
  } EmbedInfoEnum;

  /// Constructor.
  ImageStackDICOM( const Types::Coordinate tolerance = 0 /*!< Tolerance for floating point comparisons, e.g., when testing for uniform pixel/slice spacings.*/ ) : m_Tolerance( tolerance ) {}

  /// Add new DICOM image file to this stack.
  void AddImageFile( ImageFileDICOM::SmartConstPtr& image );

  /// Match new image file against this volume stack.
  bool Match ( const ImageFileDICOM& newImage /*!< New image - test whether this belongs with the ones already in this stack.*/, 
	       const Types::Coordinate numericalTolerance = 0, /*!< Numerical comparison tolerance; values with absolute difference less than this threshold are considered equal. */
	       const bool disableCheckOrientation = false /*!< Flag for disabling the checking of image orientation vectors.*/,
	       const bool ignoreAcquisitionNumber = false /*!< When this flag is set, the AcquisitionNumber DICOM tag is ignore for matching images*/ ) const;
  
  /// Write XML sidecar file.
  void WriteXML( const std::string& name /*!< Sidecar XML file name. */, const cmtk::UniformVolume& volume /*!< Previously written image volume - provides information about coordinate system etc. */,
		 const bool includeIdentifiers = false /*!< If this is set, protected "identifiers" such as device serial numbers will also be included in the XML output.*/ ) const;

  /// Write to image file.
  cmtk::UniformVolume::SmartConstPtr WriteImage ( const std::string& name /*!< File name and path for new image.*/, const Self::EmbedInfoEnum embedInfo /*!< Flag for selecting information embedded into image description.*/  ) const;

  /// Print stack information.
  void print() const;

private:
  /// Stored floating point tolerance.
  Types::Coordinate m_Tolerance;

  /// Generate custom whitespaces for XML output.
  static const char *WhitespaceWriteMiniXML( mxml_node_t*, int where);
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageStackDICOM_h_included_
