/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkGroupwiseRegistrationOutput_h_included_
#define __cmtkGroupwiseRegistrationOutput_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkGroupwiseRegistrationFunctionalBase.h>
#include <System/cmtkSmartPtr.h>
#include <Registration/cmtkReformatVolume.h>
#include <Base/cmtkInterpolator.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Class for output of groupwise registration results.
class GroupwiseRegistrationOutput
{
public:
  /// Functional base class.
  typedef GroupwiseRegistrationFunctionalBase FunctionalType;

  /// Pointer to functional base class.
  typedef FunctionalType::SmartPtr FunctionalPointer;

  /// Constructors: link to functional.
  GroupwiseRegistrationOutput( FunctionalPointer& functional = FunctionalPointer::Null ) :
    m_ExistingTemplatePath( false ),
    m_OutputRootDirectory( NULL )
  {
    this->m_Functional = functional;
  }

  /// Set flag for existing template path.
  void SetExistingTemplatePath( const bool flag )
  {
    this->m_ExistingTemplatePath = flag;
  }

  /// Set functional with implicit dynamic cast.
  template<class TFunctional>
  void SetFunctional( SmartPointer<TFunctional>& functional )
  {
    this->m_Functional = FunctionalType::SmartPtr::DynamicCastFrom( functional );
  }

  /// Set root directory for all output files.
  void SetOutputRootDirectory( const char* rootDir )
  {
    this->m_OutputRootDirectory = rootDir;
  }

  /// Write template specifications and transformations to a single file.
  bool WriteGroupwiseArchive( const char* path ) const;  
  
  /// Write each transformations to a different typedstream archive.
  bool WriteXformsSeparateArchives( const char* path, const char* templatePath );
  
  /// Reformat and write average image.
  bool WriteAverageImage( const char* path /*<! Path of output image.*/,
			  const cmtk::Interpolators::InterpolationEnum interp = cmtk::Interpolators::LINEAR /*!< Selection of interpolation method (via igsReformatVolume).*/,
			  const bool useTemplateData = false /*!< If true, template image data is included in averaging.*/ );
  
private:
  /// Pointer to functional.
  FunctionalPointer m_Functional;

  /// Flag for existing vs. generated template path.
  bool m_ExistingTemplatePath;

  /// Output root directory.
  const char* m_OutputRootDirectory;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkGroupwiseRegistrationOutput_h_included_
