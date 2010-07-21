/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkOptimizerBase_h_included_
#define __cmtkOptimizerBase_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkFunctional.h"
#include "System/cmtkCannotBeCopied.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Base class for all optimizers and meta optimizers.
class OptimizerBase :
    /// Inherit to prevent object copying.
    private CannotBeCopied
{
public:
  /// This class.
  typedef OptimizerBase Self;

  /// Return type.
  typedef Functional::ReturnType ReturnType;

  /// Parameter type.
  typedef Functional::ParameterType ParameterType;

  /// Default constructor.
  OptimizerBase() : m_FinalValue( 0.0 ) {};
  
  /// Virtual destructor.
  virtual ~OptimizerBase() {}

  /// Get final functional value.
  Self::ReturnType GetFinalValue() const
  {
    return this->m_FinalValue;
  }

protected:
  /// Set final functional value.
  void SetFinalValue( const Self::ReturnType finalValue )
  {
    this->m_FinalValue = finalValue;
  }

private:
  /// Final functional value.
  Self::ReturnType m_FinalValue;
};

//@}

} // namespace cmtk

#endif // #ifdef __cmtkOptimizerBase_h_included_
