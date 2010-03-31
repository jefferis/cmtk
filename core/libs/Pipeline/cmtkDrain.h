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

#ifndef __cmtkDrain_h_included_
#define __cmtkDrain_h_included_

#include "cmtkPipelineObject.h"

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** General data drain template class.
 */
template<class I>
class Drain : 
  public PipelineObject 
{
public:
  /// Set input object.
  void SetInput( I *const input ) 
  {   
    this->ReplaceObject( Input, input );
  }
  
  /// The actual Update() function.
  virtual long Update() 
  {
    this->CheckInputForUpdate( Input );
    return this->Superclass::Update();
  }
  
protected:
  /// Default constructor.
  Drain() { Input = NULL; }
  
  /// Destructor.
  ~Drain() { if ( Input != NULL ) Input->Delete(); }

  /// Input object.
  I *Input;

private:
  /// Convenience declaration for calls to parent class' functions.
  typedef PipelineObject Superclass;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDrain_h_included_
