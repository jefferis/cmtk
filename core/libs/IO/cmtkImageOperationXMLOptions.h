/*
//
//  Copyright 2022 SRI International
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
//  $Revision:$
//
//  $LastChangedDate:$
//
//  $LastChangedBy:$
//
*/
#ifndef __cmtkImageOperationXMLOptions_h_included_
#define __cmtkImageOperationXMLOptions_h_included_

namespace
cmtk
{
  class ImageOperationXMLOptions
  {
    public:
    const bool m_strictXML;
    const bool m_includeNDARIdentifiers;
    const bool m_includeIdentifiers;

    ImageOperationXMLOptions(
      const bool strictXML = false,
      const bool includeNDARItentifiers = false,
      const bool includeIdentifiers = false) :
        m_strictXML(strictXML),
        m_includeNDARIdentifiers(includeNDARItentifiers),
        m_includeIdentifiers(includeIdentifiers) {}
  };
}  // namespace cmtk

#endif // #ifndef __cmtkImageOperationXMLOptions_h_included_
