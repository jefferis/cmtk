/*
//
//  Copyright 2009-2010 SRI International
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

#include "cmtkSplineWarpXformITKIO.h"

#include <IO/cmtkAffineXformITKIO.h>
#include <Base/cmtkTransformChangeToSpaceAffine.h>

#include <fstream>
#include <string>
#include <typeinfo>

void
cmtk::SplineWarpXformITKIO
::Write( const std::string& filename, const SplineWarpXform& xform, const UniformVolume& refVolume, const UniformVolume& fltVolume )
{
  std::ofstream stream( filename.c_str() );
  if ( stream.good() )
    {
    // write header
    stream << "#Insight Transform File V1.0\n"
	   << "# Transform 0\n";
    
    // write ID depending on whether CMTK is using single or double precision floats for coordinates
    if ( typeid( Types::Coordinate ) == typeid( double ) )
      {
      stream << "Transform: BSplineDeformableTransform_double_3_3\n";
      }
    else
      {
      stream << "Transform: BSplineDeformableTransform_float_3_3\n";
      }

    // write parameters
    stream << "Parameters:";

    Vector3D v, vx;
    const AffineXform::SmartPtr bulkXform = xform.GetInitialAffineXform();

    for ( size_t cp = 0; cp < xform.GetNumberOfControlPoints(); ++cp )
      {
      xform.GetOriginalControlPointPositionByOffset( v, cp );
      if ( bulkXform )
	bulkXform->ApplyInPlace( v );
      xform.GetShiftedControlPointPositionByOffset( vx, cp );

      vx -= v;
      stream << " " << vx[0] << " " << vx[1] << " " << vx[2];
      }
    stream << "\n";

    // Origin of the control point grid must be transformed into physical coordinates of the reference image
    Vector3D origin( xform.m_Offset * refVolume.GetImageToPhysicalMatrix() );
    
    // Fixed parameters:
    // * Grid Size
    // * Grid Origin
    // * Grid Spacing
    // * Grid Direction
    stream << "FixedParameters: "
	   << xform.m_Dims[0] << " " << xform.m_Dims[1] << " " << xform.m_Dims[2] << " "
	   << origin[0] << " " << origin[1] << " " << origin[2] << " "
	   << xform.Spacing[0] << " " << xform.Spacing[1] << " " << xform.Spacing[2] << " "
	   << "1 0 0 0 1 0 0 0 1\n";

    if ( bulkXform )
      {
      TransformChangeToSpaceAffine toNative( *(bulkXform), refVolume, fltVolume, AnatomicalOrientationBase::SPACE_ITK );
      AffineXformITKIO::Write( stream, toNative.GetTransformation(), 1 /*idx*/ );
      }

    stream.close();
    }
}

