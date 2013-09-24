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

#include "cmtkActiveDeformationModel.h"

#include <System/cmtkConsole.h>

#include <Base/cmtkVector3D.h>
#include <Base/cmtkAffineXform.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class W>
ActiveDeformationModel<W>::ActiveDeformationModel
( const std::list< SmartPointer<W> >& deformationList, 
  const unsigned int numberOfModes,
  const bool includeScaleInModel,
  const bool includeReferenceInModel )
{
  IncludeScaleInModel = includeScaleInModel;
  IncludeReferenceInModel = includeReferenceInModel;

  unsigned int numberOfSamples = deformationList.size();
  if ( IncludeReferenceInModel )
    ++numberOfSamples;
  
  Types::Coordinate** samplePoints = Memory::ArrayC::Allocate<Types::Coordinate*>( numberOfSamples );
  unsigned int numberOfPoints = 0;
  
  typename std::list< SmartPointer<W> >::const_iterator it = deformationList.begin();

  // prepare this object to act as an actual deformation.
  this->InitGrid( (*it)->m_Domain, (*it)->m_Dims );
  // copy Origin field of first warp.
  this->m_Offset = (*it)->m_Offset;
  
  unsigned int sample = 0;
  Types::Coordinate globalScaling = 0;
  if ( IncludeReferenceInModel ) 
    {
    samplePoints[sample++] = this->MakeSamplePointsReference( *it );
    }
  
  while ( it != deformationList.end() ) 
    {
    if ( it == deformationList.begin() ) 
      {
      numberOfPoints = (*it)->m_NumberOfParameters;
      } 
    else
      {
      if ( numberOfPoints != (*it)->m_NumberOfParameters ) 
	{
	StdErr << "WARNING: differing numbers of parameters encountered in "
		  << "ActiveDeformationModel constructor. Skipping this "
		  << "sample.";
	--numberOfSamples;
	++it;
	continue;
	}
      }
    
    samplePoints[sample++] = (*it)->GetPureDeformation( this->IncludeScaleInModel );
    globalScaling += static_cast<Types::Coordinate>( log( (*it)->GetGlobalScaling() ) );
    ++it;
    }
  
  // Set Initial Affine Transform to Identity
  AffineXform::SmartPtr identity( new AffineXform() );
  this->SetInitialAffineXform( identity );
  
  // Set global scaling to average of individual scale factors, unless it
  // was preserved as part of the actual model.
  if ( ! IncludeScaleInModel ) 
    {
    this->m_GlobalScaling = exp( globalScaling / sample );
    } 
  else
    {
    this->m_GlobalScaling = 1.0;
    }

  this->Construct( samplePoints, numberOfSamples, numberOfPoints, numberOfModes );
  
  for ( unsigned int n = 0; n < numberOfSamples; ++n )
    Memory::ArrayC::Delete( samplePoints[ n ] );
  Memory::ArrayC::Delete( samplePoints );
}

template<class W>
Types::Coordinate* 
ActiveDeformationModel<W>::MakeSamplePointsReference( const W* deformation )
{
  const unsigned int numberOfParameters = deformation->m_NumberOfParameters;
  Types::Coordinate* points = Memory::ArrayC::Allocate<Types::Coordinate>( numberOfParameters );

  Types::Coordinate* ptr = points;
  for ( unsigned int pointIdx = 0; pointIdx < numberOfParameters / 3; ++pointIdx, ptr += 3 ) 
    {
    // get original (undeformed) control point position
    const Vector3D v = deformation->GetOriginalControlPointPositionByOffset( pointIdx );
    
    // copy the result into ouput array
    for ( unsigned int dim = 0; dim < 3; ++dim ) 
      ptr[dim] = v[dim];
    }
  
  return points;
}

template<class W>
Types::Coordinate* 
ActiveDeformationModel<W>::MakeSamplePoints( const W* deformation )
{
  const unsigned int numberOfParameters = deformation->m_NumberOfParameters;
  Types::Coordinate* points = Memory::ArrayC::Allocate<Types::Coordinate>( numberOfParameters );
  memcpy( points, deformation->m_Parameters, sizeof( *points ) * numberOfParameters );

  AffineXform::SmartPtr xform( deformation->GetInitialAffineXform()->MakeInverse() );
  
  if ( IncludeScaleInModel ) 
    {
    xform->SetScales( 1.0, 1.0, 1.0 );
    }
  
  Types::Coordinate* ptr = points;
  Vector3D u;
  for ( unsigned int pointIdx = 0; pointIdx < numberOfParameters / 3; ++pointIdx, ptr += 3 ) 
    {
    // undo affine transformation component
    const FixedVector<3,Types::Coordinate> v = xform->Apply( FixedVector<3,Types::Coordinate>::FromPointer( ptr ) );

    // copy the result into ouput array
    for ( unsigned int dim = 0; dim < 3; ++dim ) 
      ptr[dim] = v[dim];
    }
  
  return points;
}

template<class W>
W*
ActiveDeformationModel<W>::Compose
( const Types::Coordinate* weights )
{
  this->m_Parameters = this->Generate( this->m_Parameters, weights );

  return this;
}

template<class W>
float
ActiveDeformationModel<W>::Decompose
( const W* input, Types::Coordinate *const weights ) const
{
  CoordinateVector inputVector( this->GetNumberOfPoints(), input->GetPureDeformation( this->IncludeScaleInModel ) );
  return this->ActiveShapeModel::Decompose( &inputVector, weights );
}

/// Instance the spline warp ADM.
template class ActiveDeformationModel<SplineWarpXform>;

} // namespace cmtk
