/*
//
//  Copyright 2009 SRI International
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

#include <cmtkImageOperationCropThreshold.h>

cmtk::UniformVolume::SmartPtr
cmtk::ImageOperationCropThreshold::Apply( cmtk::UniformVolume::SmartPtr& volume )
{
  volume->AutoCrop( this->m_Threshold, true /*recrop*/ );
  
  if ( this->m_WriteRegion )
    {
    const DataGrid::RegionType& cropRegion = volume->CropRegion();
    printf( "AutoCrop %d,%d,%d,%d,%d,%d\n", cropRegion.From()[0], cropRegion.From()[1], cropRegion.From()[2], cropRegion.To()[0], cropRegion.To()[1], cropRegion.To()[2] );
    }
  
  if ( this->m_WriteXform )
    {
    cmtk::StdErr << "SORRY, this is not yet implemented!\n";
    }
  
  return cmtk::UniformVolume::SmartPtr( volume->GetCroppedVolume() );    
}
