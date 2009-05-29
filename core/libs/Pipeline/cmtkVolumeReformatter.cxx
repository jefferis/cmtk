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

#include <cmtkVolumeReformatter.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

VolumeReformatter::VolumeReformatter()
{
  ReferenceIndex = 0;
  TopIndex = 0;
  LowerThreshold = 0;
  UpperThreshold = 0;

  RescaleIndex = 0;
  RescaleOffset = 0;
  RescaleSlope = 1;

  OutputFormat = 0;
  Anonymize = 0;
  CheckerboardMode = 1;

  ImagePath = FilenamePattern = NULL;
}

VolumeReformatter::~VolumeReformatter()
{
  if ( ImagePath ) free( ImagePath );
  if ( FilenamePattern ) free( FilenamePattern );
}

void VolumeReformatter::Execute()
{
  if ( !Input[0] || !Input[1] ) return;

  UniformVolume::SmartPtr referenceVolume = Input[ReferenceIndex]->GetVolume();
  if ( ! referenceVolume ) return;

  UniformVolume::SmartPtr floatingVolume = Input[1-ReferenceIndex]->GetVolume();
  if ( ! floatingVolume ) 
    {
    return;
    }

  this->m_ReformatVolume.SetReferenceVolume( referenceVolume );
  this->m_ReformatVolume.SetFloatingVolume( floatingVolume );
  
  AffineXform::SmartPtr AffineXform = Input[1-ReferenceIndex]->GetAffineXform()->GetInverse();
  if ( ! AffineXform )
    AffineXform = Input[ReferenceIndex]->GetAffineXform();
  this->m_ReformatVolume.SetAffineXform( AffineXform );

  this->m_ReformatVolume.SetCheckerboardMode( CheckerboardMode );
  this->m_ReformatVolume.SetRescale( RescaleIndex ==ReferenceIndex, RescaleOffset, RescaleSlope );

  int midPlane = referenceVolume->GetDims()[2] / 2;
  TypedArray::SmartPtr previewData = TypedArray::SmartPtr( this->m_ReformatVolume.PlainReformat( midPlane ) );
  
  Image *output = this->GetOutput();
  output->SetDims( referenceVolume->GetDims()[0], referenceVolume->GetDims()[1] );
  output->SetSpacing( referenceVolume->GetDelta( AXIS_X, 0 ), referenceVolume->GetDelta( AXIS_Y, 0 ) );
  output->SetData( previewData );
  
  this->UpdateExecuteTime();
}

void VolumeReformatter::ExecuteReformat()
{
  this->CheckAllInputsForUpdate();
}

} // namespace cmtk
