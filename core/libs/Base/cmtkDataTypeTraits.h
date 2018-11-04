/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#ifndef __cmtkDataTypeTraits_h_included_
#define __cmtkDataTypeTraits_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkTypes.h>

#include <limits>

namespace cmtk {

/** \addtogroup Base */
//@{

/// Make signed data types.
template <class TType>
class MakeSigned {};

template <>
class MakeSigned<unsigned long int> {
 public:
  typedef signed long int Type;
};
template <>
class MakeSigned<signed long int> {
 public:
  typedef signed long int Type;
};

template <>
class MakeSigned<unsigned int> {
 public:
  typedef signed int Type;
};
template <>
class MakeSigned<signed int> {
 public:
  typedef signed int Type;
};

template <>
class MakeSigned<unsigned short int> {
 public:
  typedef signed short int Type;
};
template <>
class MakeSigned<signed short int> {
 public:
  typedef signed short int Type;
};

template <>
class MakeSigned<unsigned char> {
 public:
  typedef signed char Type;
};
template <>
class MakeSigned<signed char> {
 public:
  typedef signed char Type;
};

template <>
class MakeSigned<float> {
 public:
  typedef float Type;
};
template <>
class MakeSigned<double> {
 public:
  typedef double Type;
};

/** Data type traits */
template <class TType>
class DataTypeTraits {
 public:
  /** Scalar type ID constant for this type. */
  static const ScalarDataType DataTypeID = TYPE_NONE;
};

/** Data type traits for single-precision floating point. */
template <>
class DataTypeTraits<float> {
 public:
  /** This class. */
  typedef DataTypeTraits<float> Self;

  /// The template value type.
  typedef float ValueType;

  /** Scalar type ID constant for this type. */
  static const ScalarDataType DataTypeID = TYPE_FLOAT;

  /** Return given value converted (and rounded) to discrete type. */
  template <class T>
  static inline Self::ValueType Convert(const T value, const bool = false,
                                        const Self::ValueType = 0) {
    return static_cast<Self::ValueType>(value);
  }

  /** Return padding data value (i.e. Inf) for the given type. */
  static inline Self::ValueType ChoosePaddingValue() {
    return std::numeric_limits<Self::ValueType>::infinity();
  }

  /// Return zero value for this type.
  static inline Self::ValueType Zero() { return 0.0; }

  /// Return one value (multiplicative identity) for this type.
  static inline Self::ValueType One() { return 1.0; }
};

/** Data type traits for double-precision floating point. */
template <>
class DataTypeTraits<double> {
 public:
  /** This class. */
  typedef DataTypeTraits<double> Self;

  /// The template value type.
  typedef double ValueType;

  /** Scalar type ID constant for this type. */
  static const ScalarDataType DataTypeID = TYPE_DOUBLE;

  /** Return given value converted (and rounded) to discrete type. */
  template <class T>
  static inline Self::ValueType Convert(const T value, const bool = false,
                                        const Self::ValueType = 0) {
    return static_cast<Self::ValueType>(value);
  }

  /** Return padding data value (i.e. Inf) for the given type. */
  static inline Self::ValueType ChoosePaddingValue() {
    return std::numeric_limits<Self::ValueType>::infinity();
  }

  /// Return zero value for this type.
  static inline Self::ValueType Zero() { return 0.0; }

  /// Return one value (multiplicative identity) for this type.
  static inline Self::ValueType One() { return 1.0; }
};

/** Data type traits for unsigned char (byte). */
template <>
class DataTypeTraits<byte> {
 public:
  /** This class. */
  typedef DataTypeTraits<byte> Self;

  /// The template value type.
  typedef byte ValueType;

  /** Scalar type ID constant for this type. */
  static const ScalarDataType DataTypeID = TYPE_BYTE;

  /** Return given value converted (and rounded) to discrete type. */
  template <class T>
  static inline Self::ValueType Convert(const T value,
                                        const bool paddingFlag = false,
                                        const Self::ValueType paddingData = 0) {
    using namespace std;
    if (MathUtil::IsFinite(value)) {
      return (Self::ValueType)(
          (value < std::numeric_limits<Self::ValueType>::min())
              ? std::numeric_limits<Self::ValueType>::min()
              : (value + 0.5 > std::numeric_limits<Self::ValueType>::max())
                    ? std::numeric_limits<Self::ValueType>::max()
                    : floor(value + 0.5));
    } else {
      if (paddingFlag)
        return paddingData;
      else
        return ChoosePaddingValue();
    }
  }

  /** Return padding data value for the given type. */
  static inline Self::ValueType ChoosePaddingValue() { return 255; }

  /// Return zero value for this type.
  static inline Self::ValueType Zero() { return 0; }

  /// Return one value (multiplicative identity) for this type.
  static inline Self::ValueType One() { return 1; }
};

/** Data type traits for signed char. */
template <>
class DataTypeTraits<char> {
 public:
  /** This class. */
  typedef DataTypeTraits<char> Self;

  /// The template value type.
  typedef char ValueType;

  /** Scalar type ID constant for this type. */
  static const ScalarDataType DataTypeID = TYPE_CHAR;

  /** Return given value converted (and rounded) to discrete type. */
  template <class T>
  static inline Self::ValueType Convert(const T value,
                                        const bool paddingFlag = false,
                                        const Self::ValueType paddingData = 0) {
    using namespace std;
    if (MathUtil::IsFinite(value)) {
      return (Self::ValueType)(
          (value < std::numeric_limits<Self::ValueType>::min())
              ? std::numeric_limits<Self::ValueType>::min()
              : (value + 0.5 > std::numeric_limits<Self::ValueType>::max())
                    ? std::numeric_limits<Self::ValueType>::max()
                    : floor(value + 0.5));
    } else {
      if (paddingFlag)
        return paddingData;
      else
        return ChoosePaddingValue();
    }
  }

  /** Return padding data value for the given type. */
  static inline Self::ValueType ChoosePaddingValue() { return -1; }

  /// Return zero value for this type.
  static inline Self::ValueType Zero() { return 0; }

  /// Return one value (multiplicative identity) for this type.
  static inline Self::ValueType One() { return 1; }
};

/** Data type traits for signed short. */
template <>
class DataTypeTraits<signed short> {
 public:
  /** This class. */
  typedef DataTypeTraits<signed short> Self;

  /// The template value type.
  typedef signed short ValueType;

  /** Scalar type ID constant for this type. */
  static const ScalarDataType DataTypeID = TYPE_SHORT;

  /** Return given value converted (and rounded) to discrete type. */
  template <class T>
  static inline Self::ValueType Convert(const T value,
                                        const bool paddingFlag = false,
                                        const Self::ValueType paddingData = 0) {
    using namespace std;
    if (MathUtil::IsFinite(value)) {
      return (Self::ValueType)(
          ((value < std::numeric_limits<Self::ValueType>::min())
               ? std::numeric_limits<Self::ValueType>::min()
               : (value + 0.5 > std::numeric_limits<Self::ValueType>::max())
                     ? std::numeric_limits<Self::ValueType>::max()
                     : floor(value + 0.5)));
    } else {
      if (paddingFlag)
        return paddingData;
      else
        return ChoosePaddingValue();
    }
  }

  /** Return padding data value for the given type. */
  static inline Self::ValueType ChoosePaddingValue() { return -1; }

  /// Return zero value for this type.
  static inline Self::ValueType Zero() { return 0; }

  /// Return one value (multiplicative identity) for this type.
  static inline Self::ValueType One() { return 1; }
};

/** Data type traits for unsigned short. */
template <>
class DataTypeTraits<unsigned short> {
 public:
  /** This class. */
  typedef DataTypeTraits<unsigned short> Self;

  /// The template value type.
  typedef unsigned short ValueType;

  /** Scalar type ID constant for this type. */
  static const ScalarDataType DataTypeID = TYPE_USHORT;

  /** Return given value converted (and rounded) to discrete type. */
  template <class T>
  static inline Self::ValueType Convert(const T value,
                                        const bool paddingFlag = false,
                                        const Self::ValueType paddingData = 0) {
    using namespace std;
    if (MathUtil::IsFinite(value)) {
      return (Self::ValueType)(
          (value < std::numeric_limits<Self::ValueType>::min())
              ? std::numeric_limits<Self::ValueType>::min()
              : (value + 0.5 > std::numeric_limits<Self::ValueType>::max())
                    ? std::numeric_limits<Self::ValueType>::max()
                    : floor(value + 0.5));
    } else {
      if (paddingFlag)
        return paddingData;
      else
        return ChoosePaddingValue();
    }
  }

  /** Return padding data value for the given type. */
  static inline Self::ValueType ChoosePaddingValue() {
    return std::numeric_limits<Self::ValueType>::max();
  }

  /// Return zero value for this type.
  static inline Self::ValueType Zero() { return 0; }

  /// Return one value (multiplicative identity) for this type.
  static inline Self::ValueType One() { return 1; }
};

/** Data type traits for int. */
template <>
class DataTypeTraits<int> {
 public:
  /** This class. */
  typedef DataTypeTraits<int> Self;

  /// The template value type.
  typedef int ValueType;

  /** Scalar type ID constant for this type. */
  static const ScalarDataType DataTypeID = TYPE_INT;

  /** Return given value converted (and rounded) to discrete type. */
  template <class T>
  static inline Self::ValueType Convert(const T value,
                                        const bool paddingFlag = false,
                                        const Self::ValueType paddingData = 0) {
    using namespace std;
    if (MathUtil::IsFinite(value)) {
      return (Self::ValueType)(
          (value < std::numeric_limits<Self::ValueType>::min())
              ? std::numeric_limits<Self::ValueType>::min()
              : (value + 0.5 > std::numeric_limits<Self::ValueType>::max())
                    ? std::numeric_limits<Self::ValueType>::max()
                    : floor(value + 0.5));
    } else {
      if (paddingFlag)
        return paddingData;
      else
        return ChoosePaddingValue();
    }
  }

  /** Return padding data value for the given type. */
  static inline Self::ValueType ChoosePaddingValue() { return -1; }

  /// Return zero value for this type.
  static inline Self::ValueType Zero() { return 0; }

  /// Return one value (multiplicative identity) for this type.
  static inline Self::ValueType One() { return 1; }
};

/** Data type traits for unsigned int. */
template <>
class DataTypeTraits<unsigned int> {
 public:
  /** This class. */
  typedef DataTypeTraits<unsigned int> Self;

  /// The template value type.
  typedef unsigned int ValueType;

  /** Scalar type ID constant for this type. */
  static const ScalarDataType DataTypeID = TYPE_UINT;

  /** Return given value converted (and rounded) to discrete type. */
  template <class T>
  static inline Self::ValueType Convert(const T value,
                                        const bool paddingFlag = false,
                                        const Self::ValueType paddingData = 0) {
    using namespace std;
    if (MathUtil::IsFinite(value)) {
      return (Self::ValueType)(
          (static_cast<Self::ValueType>(value) <
           std::numeric_limits<Self::ValueType>::min())
              ? std::numeric_limits<Self::ValueType>::min()
              : (value + 0.5 > std::numeric_limits<Self::ValueType>::max())
                    ? std::numeric_limits<Self::ValueType>::max()
                    : floor(value + 0.5));
    } else {
      if (paddingFlag)
        return paddingData;
      else
        return ChoosePaddingValue();
    }
  }

  /** Return padding data value for the given type. */
  static inline Self::ValueType ChoosePaddingValue() {
    return static_cast<Self::ValueType>(-1);
  }

  /// Return zero value for this type.
  static inline Self::ValueType Zero() { return 0; }

  /// Return one value (multiplicative identity) for this type.
  static inline Self::ValueType One() { return 1; }
};

//@}

}  // namespace cmtk

#endif  // #ifndef __cmtkDataTypeTraits_h_included_
