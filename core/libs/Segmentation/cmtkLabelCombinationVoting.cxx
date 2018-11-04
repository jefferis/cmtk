/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011, 2014 SRI International
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

#include <Segmentation/cmtkLabelCombinationVoting.h>

#include <algorithm>
#include <vector>

namespace cmtk {

/** \addtogroup Segmentation */
//@{

LabelCombinationVoting::LabelCombinationVoting(
    const std::vector<TypedArray::SmartPtr> &data) {
  const size_t nValues = data[0]->GetDataSize();
  this->m_Result =
      TypedArray::SmartPtr(TypedArray::Create(TYPE_SHORT, nValues));
  this->m_Result->SetDataClass(DATACLASS_LABEL);

  size_t numberOfClasses = 1;
  for (size_t k = 0; k < data.size(); ++k) {
    const Types::DataItemRange range = data[k]->GetRange();
    numberOfClasses =
        std::max(numberOfClasses, 1 + static_cast<size_t>(range.m_UpperBound));
  }

  std::vector<unsigned int> label(1 + numberOfClasses);

  for (size_t i = 0; i < nValues; ++i) {
    std::fill(label.begin(), label.end(), 0);

    for (size_t curr = 0; curr < data.size(); ++curr) {
      Types::DataItem v;
      if (data[curr]->Get(v, i)) {
        ++label[std::min(numberOfClasses, static_cast<size_t>(v))];
      }
    }

    // Compute winner of label voting.

    short maxLab = 0;
    unsigned int maxCnt = 0;

    for (size_t lab = 0; lab < numberOfClasses; ++lab) {
      // do something with tie case
      if (label[lab] > maxCnt) {
        maxCnt = label[lab];
        maxLab = lab;
      } else {
        if (label[lab] == maxCnt) {
          maxLab = -1;
        }
      }
    }

    this->m_Result->Set(maxLab, i);
  }
}

}  // namespace cmtk
