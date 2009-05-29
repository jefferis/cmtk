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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkLabelCombinationVoting.h>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

LabelCombinationVoting::LabelCombinationVoting( const std::vector<TypedArray::SmartPtr>& data )
{
  const size_t nValues = data[ 0 ]->GetDataSize();
  m_Result = TypedArray::SmartPtr( TypedArray::Create( TYPE_SHORT, nValues ) );
  
#pragma omp parallel for  
  for ( size_t i = 0; i < nValues; ++i )
    {
    short label[256];
    memset( label, 0, sizeof( label ) );

    for ( size_t curr = 0; curr < data.size(); ++curr )
      {
      Types::DataItem v;
      if ( data[ curr ]->Get( v, i ) ) 
        {
        ++label[ static_cast<byte>( v ) ];
        }
      }

    // Compute winner of label voting.

    int maxLab = 0;
    int maxCnt = 0;
   
    for ( int lab=0; lab < 256; ++lab ) 
      {
      // do something with tie case
      if ( label[ lab ] > maxCnt ) 
        {
        maxCnt = label[ lab ];
        maxLab = lab;
        } 
      else
	{
	if ( label[lab] == maxCnt )
	  {
	  maxLab = -1;
	  }
	}
      }
  
    this->m_Result->Set( maxLab, i ); 
    }
}

} // namespace cmtk
