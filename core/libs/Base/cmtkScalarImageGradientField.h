/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#ifndef __cmtkScalarImageGradientField_h_included_
#define __cmtkScalarImageGradientField_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkImageTemplate.h>
#include <Base/cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Compute gradient field from scalar image
class ScalarImageGradientField
{
public:
  /// This class.
  typedef ScalarImageGradientField Self;

  /// The gradient field type.
  typedef ImageTemplate< FixedVector<3,Types::Coordinate> > GradientFieldType;

  /// Constructor: compute gradient field
  ScalarImageGradientField( const UniformVolume& volume /*!< Scalar image volume to compute gradients from. */ );

  /// Get gradient field.
  Self::GradientFieldType::SmartPtr Get()
  {
    return this->m_GradientField;
  }

private:
  /// The gradient field.
  Self::GradientFieldType::SmartPtr m_GradientField;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkScalarImageGradientField_h_included_
