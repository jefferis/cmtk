/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include "cmtkLabelCombinationLocalWeighting.h"

#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkRegionIndexIterator.h>

#include <Registration/cmtkTypedArraySimilarity.h>

#ifdef _OPENMP
#  include <omp.h>
#endif

void
cmtk::LabelCombinationLocalWeighting::AddAtlasImage
( const UniformVolume::SmartConstPtr image )
{
  if ( !this->m_TargetImage->GridMatches( *image ) )
    {
    StdErr << "Atlas intensity image grid does not match target image.\n";
    throw ExitException( 1 );
    }

  this->m_AtlasImages.push_back( image );
}