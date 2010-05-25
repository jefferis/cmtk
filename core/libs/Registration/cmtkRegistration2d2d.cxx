/*
//
//  Copyright 2004-2010 SRI International
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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkRegistration2d2d.h>

#include <cmtkFunctionalAffine2D.h>
#include <cmtkOptimizer.h>
#include <cmtkBestNeighbourOptimizer.h>

#include <cmtkVector.h>
#include <cmtkTypedArraySimilarity.h>

#include <cmtkProtocolCallback.h>

#include <cmtkPGM.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

void 
Registration2d2d::Register
( CoordinateMatrix3x3& matrix, ScalarImage::SmartPtr& refImage,
  ScalarImage::SmartPtr& fltImage, const ScalarImage::RegionType* fltROI )
{
  ScalarImage::SmartPtr roi( new ScalarImage( *(fltImage) ) );

  if ( fltROI ) 
    {
    roi->SetROI( *fltROI );
    Types::Coordinate v[8];
    matrix.Decompose( v );
    v[0] += fltROI->From()[AXIS_X] * roi->GetPixelSize( AXIS_X );
    v[1] += fltROI->From()[AXIS_Y] * roi->GetPixelSize( AXIS_Y );
    matrix.Compose( v );
    }

  Register( matrix, refImage, roi );
}

void 
Registration2d2d::Register
( CoordinateMatrix3x3& matrix, ScalarImage::SmartPtr& refImage,
  ScalarImage::SmartPtr& fltImage )
{
  SmartPointer<FunctionalAffine2D> functional
    ( new FunctionalAffine2D( refImage, fltImage ) );
  functional->SetSimilarityMeasure( ScalarImageSimilarity::MI );
    
  CoordinateMatrix3x3 init( matrix );
  
  BestNeighbourOptimizer optimizer;
  //  RegistrationCallback_P callback
  //    ( new ProtocolCallback( "track.txt", true ) );
  //  optimizer.SetCallback( callback );
    
  optimizer.SetFunctional( Functional::SmartPtr::DynamicCastFrom( functional ) );
  CoordinateVector v( 8 );
  matrix.Decompose( v.Elements );

  const float explore1 = 5.0;
  const float accuracy1 = 0.5;
  functional->SetNumberDOFs( 3 );
  optimizer.Optimize( v, explore1, accuracy1 );

  //  const float explore2 = 2.0;
  //  const float accuracy2 = 0.125;
  //  functional->SetNumberDOFs( 3 );
  //  optimizer.Optimize( v, explore1, accuracy2 );

  matrix.Compose( v.Elements );

  ScalarImage::SmartPtr finalImage( refImage->InterpolateFrom( fltImage, &matrix ) );

  static int cnt = 0;
  char fname[12];
  snprintf( fname, sizeof( fname ), "fin%03d.pgm", cnt++ );
  PGM::Write16bit( fname, finalImage );
}

} // namespace cmtk
