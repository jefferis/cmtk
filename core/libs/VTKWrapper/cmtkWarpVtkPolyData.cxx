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

#include <cmtkWarpVtkPolyData.h>

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

namespace
cmtk
{

/** \addtogroup VTKWrapper */
//@{

void
WarpVtkPolyData::SetWarpXform( WarpXform::SmartPtr& warpXform )
{
  if ( this->m_WarpXform == warpXform ) return;
  this->m_WarpXform = warpXform;
  this->Modified();
}

void
WarpVtkPolyData::SetAffineXform( AffineXform::SmartPtr& affineXform )
{
  if ( this->m_AffineXform == affineXform ) return;
  this->m_AffineXform = affineXform;
  this->Modified();
}

void WarpVtkPolyData::Execute()
{
  vtkPolyData *input = this->GetInput();
  vtkPolyData *output = this->GetOutput();

  vtkPoints *inPoints = input->GetPoints();
  vtkDataArray *inCoordArray = inPoints->GetData();

  if ( !this->m_WarpXform && !this->m_AffineXform ) 
    {
    output->SetPoints( inPoints );
    } 
  else
    {
    float *inPointPtr = dynamic_cast<vtkFloatArray*>( inCoordArray )->GetPointer( 0 );
    
    int numVectors = inCoordArray->GetNumberOfTuples();
    
    vtkDataArray *outCoordArray = vtkFloatArray::New();
    outCoordArray->SetNumberOfComponents( 3 );
    float *outPointPtr = dynamic_cast<vtkFloatArray*>( outCoordArray )->WritePointer( 0, 3*numVectors );
    
    Vector3D v;
    if ( this->m_WarpXform ) 
      {
      for ( int idx = 0; idx<numVectors; ++idx ) 
	{
	v = Vector3D( inPointPtr+3*idx );
	this->m_WarpXform->ApplyInPlace( v );
	outPointPtr[3*idx] = v[0];
	outPointPtr[3*idx+1] = v[1];
	outPointPtr[3*idx+2] = v[2];
	}
      } 
    else
      {
      for ( int idx = 0; idx<numVectors; ++idx ) 
	{
	v = Vector3D( inPointPtr+3*idx );
	this->m_AffineXform->ApplyInPlace( v );
	outPointPtr[3*idx] = v[0];
	outPointPtr[3*idx+1] = v[1];
	outPointPtr[3*idx+2] = v[2];
	}
      }
    
    vtkPoints *outPoints = vtkPoints::New();
    outPoints->SetData( outCoordArray );
    outCoordArray->Delete();
    output->SetPoints( outPoints );
    outPoints->Delete();
    }
  
  vtkCellArray *cellArray = input->GetPolys();
  output->SetPolys( cellArray );

  // transfer point data, i.e., scalars, normals, gradients.
  vtkPointData *inPointData = input->GetPointData();
  vtkPointData *outPointData = output->GetPointData();
  if ( inPointData ) 
    {
    outPointData->SetScalars( inPointData->GetScalars() );
    outPointData->SetVectors( inPointData->GetVectors() );
    outPointData->SetNormals( inPointData->GetNormals() );
    }
}

} // namespace cmtk
