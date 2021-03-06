/*
//
//  Copyright 2016 Google, Inc.
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include "cmtkGroupwiseRegistrationFunctionalIO.h"

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkSplineWarpXform.h>
#include <IO/cmtkClassStreamAffineXform.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

ClassStreamOutput&
operator<<
  ( ClassStreamOutput& stream, const GroupwiseRegistrationFunctionalBase& func )
{
  const UniformVolume* templateGrid = func.GetTemplateGrid();
  stream.Begin( "template" );
  stream.WriteIntArray( "dims", templateGrid->GetDims().begin(), 3 );
  stream.WriteCoordinateArray( "delta", templateGrid->Deltas().begin(), 3 );
  stream.WriteCoordinateArray( "size", templateGrid->m_Size.begin(), 3 );
  stream.WriteCoordinateArray( "origin", templateGrid->m_Offset.begin(), 3 );
  stream.End();
  
  for ( size_t idx = 0; idx < func.GetNumberOfTargetImages(); ++idx )
    {
    const UniformVolume* target = func.GetOriginalTargetImage( idx );
    stream.WriteString( "target", target->GetMetaInfo( META_FS_PATH ).c_str() );
    
    const Xform* xform = func.GetGenericXformByIndex( idx );
    
    const AffineXform* affineXform = dynamic_cast<const AffineXform*>( xform );
    if ( affineXform )
      stream << (*affineXform);

    const SplineWarpXform* splineXform = dynamic_cast<const SplineWarpXform*>( xform );
    if ( splineXform )
      stream << splineXform;
    }
  
  return stream;
}

} // namespace cmtk
