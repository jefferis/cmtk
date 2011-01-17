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

#ifndef __cmtkFilter_h_included_
#define __cmtkFilter_h_included_

#include <Pipeline/cmtkSource.h>
#include <Pipeline/cmtkPipelineObject.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Basic filter class.
 * This class combines the data source functions inherited from Source
 * with an additional input port. It therefore serves as a template for
 * all classes transforming an input into an output object. Both, input and 
 * output type are defined by template parameters "I" and "O", respectively.
 * "O" is passed directly to the Source parent class.
 *\see Source
 */
template<class I, class O> class Filter : public Source<O> 
{
public:
  /// Replace the current Input object with a new one.
  void SetInput ( I *const input ) 
  {
    this->ReplaceObject( Input, input );
  }

  /** Update this object.
   * Check for changes in the Input object first, then call inherited Update()
   * function from Object.
   *\see Object#Update
   */
  virtual long Update () 
  {
    this->CheckInputForUpdate( Input );
    return this->PipelineObject::Update();
  }

protected:
  /// Default constructor.
  Filter() { Input = NULL; }

  /** Destructor.
   * Unregister from the Input object if one was set.
   */
  virtual ~Filter() 
  { 
    if ( Input ) 
      Input->Delete();
  }

  /// The actual input object.
  I *Input;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFilter_h_included_
