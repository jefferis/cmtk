/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#include <cmtkDataGridMorphologicalOperators.h>

#include <map>
#include <vector>

cmtk::TypedArray::SmartPtr
cmtk::DataGridMorphologicalOperators::GetConnectedComponents( const bool sortBySize ) const
{
  const size_t numberOfPixels = this->m_DataGrid->GetNumberOfPixels();

  // this vector will hold the component index per pixel; since we move strictly forward and only look back in the algorithm below,
  // there is no need to initialize the vector elements.
  std::vector<int> result( numberOfPixels);

  const DataGrid::IndexType& dims = this->m_DataGrid->GetDims();

  DataGrid::IndexType relative;
  relative[0] = this->m_DataGrid->GetNextI();
  relative[1] = this->m_DataGrid->GetNextJ();
  relative[2] = this->m_DataGrid->GetNextK();

  std::map<int,int> linkMap;
  int nextComponent = 0;

  DataGrid::IndexType index;
  size_t offset = 0;
  for ( index[2] = 0; index[2] < dims[2]; ++index[2] )
    {
    for ( index[1] = 0; index[1] < dims[1]; ++index[1] )
      {
      for ( index[0] = 0; index[0] < dims[0]; ++index[0], ++offset )
	{
	// initialize as "background"
	int component = 0;
	  
	// foreground pixel?
	if ( this->m_DataGrid->GetDataAt( offset ) != 0 )
	  {
	  // loop over x,y,z neighbor
	  for ( int dim = 0; dim < 3; ++dim )
	    {
	    // is there a preceding neighbor in this direction?
	    if ( index[dim] )
	      {
	      // get component ID for that neighbor
	      const int existing = result[offset - relative[dim]];
	      // is there something there?
	      if ( existing )
		{
		// did we already get a different component ID from another neighbor?
		if ( component > existing ) // note: this implies "component != 0" because "existing" is at least 0.
		  {
		  // link old and new component via this pixel
		  linkMap[component] = existing;
		  }
		// mark current pixel as belonging to the same component
		component = existing;
		}
	      }
	    }
	  
	  // none of the neighbors is foreground, so use next available component ID.
	  if ( !component )
	    {
	    component = ++nextComponent;
	    }
	  }
	
	// set component result for this pixel
	result[offset] = component;
	}
      }
    }
  
  // now collapse all component indexes into their unique representative
  for ( int component = 0; component < nextComponent; ++component )
    {
    int mapTo = component;
    while ( (linkMap.find( mapTo ) != linkMap.end()) && (mapTo != linkMap[mapTo]) )
      mapTo = linkMap[mapTo];
    
    linkMap[component] = mapTo;
    }
  
  // re-number components
  TypedArray::SmartPtr resultArray( TypedArray::Create( TYPE_INT, numberOfPixels ) );
  for ( size_t px = 0; px < numberOfPixels; ++px )
    {
    resultArray->Set( linkMap[result[px]], px );
    }
  
  return resultArray;
}
