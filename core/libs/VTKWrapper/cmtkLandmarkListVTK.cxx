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

#include <cmtkLandmarkListVTKVTK.h>

namespace
cmtk
{

/** \addtogroup VTKWrapper */
//@{

vtkPoints*
LandmarkListVTK::GetVtkPoints() const
{
  vtkPoints* points = vtkPoints::New();

  const_iterator it = this->begin();
  while ( it != this->end() ) 
    {
    points->InsertNextPoint( (*it)->GetLocation() );
    ++it;
    }
  
  return points;
}

vtkPoints* 
LandmarkListVTK::GetMatchedVtkPoints( vtkPoints*& targetPoints, const LandmarkListVTK *targetLL ) const
{
  vtkPoints* points = vtkPoints::New();
  targetPoints = vtkPoints::New();
  
  const_iterator it = this->begin();
  while ( it != this->end() ) 
    {
    const SmartPointer<Landmark> targetLM = targetLL->FindByName( (*it)->GetName() );
    if ( targetLM.GetPtr() )
      {
      points->InsertNextPoint( (*it)->GetLocation() );
      targetPoints->InsertNextPoint( targetLM->GetLocation() );
      }
    ++it;
    }
  
  return points;
}

} // namespace cmtk
