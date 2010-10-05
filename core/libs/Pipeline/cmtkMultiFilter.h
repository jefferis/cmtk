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

#ifndef __cmtkMultiFilter_h_included_
#define __cmtkMultiFilter_h_included_

#include <cmtkconfig.h>

#include <Pipeline/cmtkSource.h>
#include <Pipeline/cmtkPipelineObject.h>

#include <list>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Filter with several inputs.
 * This class combines the data source functions inherited from Source
 * with update control for an arbitrary number of input port. It therefore 
 * serves as a template for all classes transforming more than one input into
 * an output object. For just one input, Filter is probably more efficient
 * as it gets along without the STL "list" class.
 *@see Source
 *@see Filter
 */
template<class O> 
class MultiFilter : 
    public Source<O> 
{
public:
  template<class I> void RegisterInput( I** input ) 
  {
    if ( input ) 
      {
      this->m_InputList.push_back( (PipelineObject**) input );
      }
  }
  
  template<class I> void UnregisterInput( const I** input ) 
  {
    if ( input ) 
      {
      InputListType::iterator it = this->m_InputList.begin();
      while ( it != this->m_InputList.end() ) 
	{
	if ( *it == input ) 
	  {
	  this->m_InputList.erase( it );
	  }
	++it;
	}
      }
  }
  
  /** Update this object.
   * Check for changes in all input objects first, then call inherited Update()
   * function from PipelineObject.
   *@see PipelineObject#Update
   */
  virtual long Update () 
  {
    InputListType::iterator it = this->m_InputList.begin();
    while ( it != this->m_InputList.end() ) 
      {
      if ( **it ) 
	this->CheckInputForUpdate( **it );
      ++it;
      }
    return this->PipelineObject::Update();
  }

protected:
  /// Default constructor.
  MultiFilter() {}

  /** Destructor.
   * Empty list of input objects.
   */
  virtual ~MultiFilter() 
  { 
    while ( ! this->m_InputList.empty() )
      this->m_InputList.pop_back();
  }
  
  /// Type for the STL list holding pointers to PipelineObjects.
  typedef std::list<PipelineObject**> InputListType;

  /// The actual input object.
  InputListType m_InputList;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMultiFilter_h_included_
