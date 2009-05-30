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

#ifndef __cmtkSource_h_included_
#define __cmtkSource_h_included_

#include <cmtkconfig.h>

#include <cmtkPipelineObject.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{
/** General data source object.
 * This class provides common fields and acces functions for data source
 * objects. The output object type is defined by the template parameter "O".
 */
template<class O>
class Source : 
  /// This is a pipeline object.
  public PipelineObject 
{
protected:
  /// Default constructor.
  Source() { Output = NULL; }

  /** Destructor.
   * Unregister from the output object.
   */
  virtual ~Source() 
  { 
    if ( Output ) Output->Unregister( this );
  }
  
public:
  /// Get virtual class name.
  virtual const char *GetClassName() const 
  { 
    return "Source"; 
  }

  /** Get output object.
   * If the output object does not yet exist, a new objet is created. The
   * current Source object is then registered as the new output object's
   * primary owner.
   */
  virtual O *GetOutput()
  {
    if ( Output == NULL ) 
      {
      Output = O::New();
      Output->Register( this );
      }
    return Output;
  }
  
protected:
  /// Pointer to the actual output object.
  O *Output;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSource_h_included_
