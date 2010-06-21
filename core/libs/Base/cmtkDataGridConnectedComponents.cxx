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

#include <cmtkUnionFind.h>

#include <map>
#include <vector>
#include <list>

cmtk::TypedArray::SmartPtr
cmtk::DataGridMorphologicalOperators::GetBinaryConnectedComponents() const
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

  UnionFind<int> connected;
  int nextComponent = 1;

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
	  for ( int dim = 2; dim >= 0; --dim )
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
		if ( component && (component != existing) ) // note: this implies "component != 0" because "existing" is at least 1 in this branch.
		  {
		  // link old and new component via this pixel
		  connected.Union( connected.Find( component ), connected.Find( existing ) );
		  }
		// mark current pixel as belonging to the same component
		component = existing;
		}
	      }
	    }
	  
	  // none of the neighbors is foreground, so use next available component ID.
	  if ( !component )
	    {
	    component = nextComponent++;
	    connected.Insert( component );
	    }
	  }
	
	// set component result for this pixel
	result[offset] = component;
	}
      }
    }
  
  // now collapse all component indexes into their unique representative
  std::map<int,int> linkMap;
  for ( int component = 1; component < nextComponent; ++component )
    {
    linkMap[component] = connected.FindKey( component );
    }
  
  // re-number components
  TypedArray::SmartPtr resultArray( TypedArray::Create( TYPE_INT, numberOfPixels ) );
#pragma omp parallel for
  for ( size_t px = 0; px < numberOfPixels; ++px )
    {
    resultArray->Set( linkMap[result[px]], px );
    }
  
  return resultArray;
}

cmtk::TypedArray::SmartPtr 
cmtk::DataGridMorphologicalOperators::GetRegionsRenumberedBySize() const
{
  const size_t numberOfPixels = this->m_DataGrid->GetNumberOfPixels();
  
  // first, count pixels in each region
  std::map<int,int> regionSizeMap;
  for ( size_t px = 0; px < numberOfPixels; ++px )
    {
    const int value = static_cast<int>( this->m_DataGrid->GetDataAt( px ) );
    if ( value )
      ++regionSizeMap[value];
    }

  // now create list sorted by region size
  std::list< std::pair<const int,int> > sortedList;
  for ( std::map<int,int>::const_iterator it = regionSizeMap.begin(); it != regionSizeMap.end(); ++it )
    {
    std::list< std::pair<const int,int> >::iterator ins = sortedList.begin();
    while ( (ins != sortedList.end()) && (ins->second >= it->second) )
      ++ins;
    
    sortedList.insert( ins, *it ); // need to explicitly create new pair because some STL implementation have it->first as "const".
    }

  // create renumbering lookup map
  std::map<int,int> renumberMap;
  int component = 1;
  for ( std::list< std::pair<const int,int> >::iterator it = sortedList.begin(); it != sortedList.end(); ++it )
    {
    renumberMap[it->first] = component++;
    }
  
  // re-number components
  TypedArray::SmartPtr resultArray( TypedArray::Create( TYPE_INT, numberOfPixels ) );
  for ( size_t px = 0; px < numberOfPixels; ++px )
    {
    resultArray->Set( renumberMap[static_cast<int>( this->m_DataGrid->GetDataAt( px ) )], px );
    }
  
  return resultArray;  
}
