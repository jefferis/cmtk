/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkUnits_h_included_
#define __cmtkUnits_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

namespace
Units
{

/// Angle of rotation in degrees.
class UnitBase
{
public:
  /// Constructor.
  explicit UnitBase( const double value ) : m_Value( value ) {};

  /// Get value.
  double Value() const
  {
    return this->m_Value;
  }

  /// Constant: pi.
  static double Pi()
  {
    return 3.14159265358979323846;
  }
  
private:
  /// Actual value.
  double m_Value;
};

/// Template for arithmetic on units (treated as a vector space)
template<class T>
class Arithmetic
{
public:
  /// Left-hand scalar multiplication.
  friend const T operator*( const double lhs, const T& rhs )
  {
    return T( lhs * rhs.Value() );
  }

  /// Right-hand scalar multiplication.
  friend const T operator*( const T& lhs, const double rhs )
  {
    return T( rhs * lhs.Value() );
  }

  /// Addition.
  friend const T operator*( const T& lhs, const T& rhs )
  {
    return T( lhs.Value() + rhs.Value() );
  }
};

/// Forward declaration.
class Radians;

/// Angle of rotation in degrees.
class Degrees :
    public UnitBase, public Arithmetic<Degrees>
{
public:
  /// Constructor.
  explicit Degrees( const double value = 0 ) : UnitBase( value ) {};

  /// Conversion constructor.
  inline Degrees( const Radians& radians );
};

/// Angle of rotation in radians.
class Radians :
    public UnitBase, public Arithmetic<Radians>
{
public:
  /// Constructor.
  explicit Radians( const double value = 0 ) : UnitBase( value ) {};

  /// Conversion constructor.
  inline Radians( const Degrees& degrees );
};

inline Degrees::Degrees( const Radians& radians ) : UnitBase( radians.Value() / UnitBase::Pi() * 180 ) {};
inline Radians::Radians( const Degrees& degrees ) : UnitBase( degrees.Value() * UnitBase::Pi() / 180 ) {};

/// Forward declaration.
class GaussianFWHM;

/// Parameter "\sigma" of Gaussian kernel
class GaussianSigma :
    public UnitBase, public Arithmetic<GaussianSigma>
{
public:
  /// Constructor.
  explicit GaussianSigma( const double value = 0 ) : UnitBase( value ) {};

  /// Conversion constructor.
  inline GaussianSigma( const GaussianFWHM& radians );
};

/// Full width at half maximum of Gaussian kernel.
class GaussianFWHM :
    public UnitBase, public Arithmetic<GaussianFWHM>
{
public:
  /// Constructor.
  explicit GaussianFWHM( const double value = 0 ) : UnitBase( value ) {};

  /// Conversion constructor.
  inline GaussianFWHM( const GaussianSigma& degrees );
};

inline GaussianSigma::GaussianSigma( const GaussianFWHM& fwhm ) : UnitBase( fwhm.Value() / 2.354820045 ) {};
inline GaussianFWHM::GaussianFWHM( const GaussianSigma& sigma ) : UnitBase( sigma.Value() * 2.354820045 ) {};

}

//@}

} // namespace cmtk

#endif // #define __cmtkUnits_h_included_
