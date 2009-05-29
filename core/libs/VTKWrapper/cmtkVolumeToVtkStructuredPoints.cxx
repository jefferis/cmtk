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

#include <cmtkVolumeToVtkStructuredPoints.h>

#include <vtkStructuredPoints.h>
#include <vtkPointData.h>
#include <vtkCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkShortArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkObjectFactory.h>

#include <cmtkTypedArray.h>

namespace
cmtk
{

/** \addtogroup VTKWrapper */
//@{

VolumeToVtkStructuredPoints::VolumeToVtkStructuredPoints () :
  m_Volume( NULL ),
  m_AffineXform( NULL )
{}

VolumeToVtkStructuredPoints::VolumeToVtkStructuredPoints
(const VolumeToVtkStructuredPoints& other)
{
  this->m_Volume = other.m_Volume;
  this->m_AffineXform = other.m_AffineXform;
}

void
VolumeToVtkStructuredPoints::operator=
(const VolumeToVtkStructuredPoints& other) 
{
  this->m_Volume = other.m_Volume;
  this->m_AffineXform = other.m_AffineXform;
}

VolumeToVtkStructuredPoints* VolumeToVtkStructuredPoints::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("VolumeToVtkStructuredPoints");
  if (ret) 
    {
    return (VolumeToVtkStructuredPoints*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new VolumeToVtkStructuredPoints;
}

void
VolumeToVtkStructuredPoints::SetVolume( UniformVolume::SmartPtr& volume )
{
  this->Modified();
  if ( this->m_Volume == volume ) return;
  this->m_Volume = volume;
}

void
VolumeToVtkStructuredPoints::SetXform
( AffineXform::SmartPtr& affineXform )
{
  this->Modified();
  if ( this->m_AffineXform == affineXform ) return;
  this->m_AffineXform = affineXform;
}

void VolumeToVtkStructuredPoints::Execute()
{
  if ( ! this->m_Volume ) return;

  TypedArray::SmartPtr volumeData = this->m_Volume->GetData();
  if ( ! volumeData ) return;
  
  vtkStructuredPoints *output = this->GetOutput();
  
  output->SetDimensions( this->m_Volume->GetDims() );

  double ar[3];
  for ( int idx=0; idx<3; ++idx )
    ar[idx] = this->m_Volume->GetDelta( idx, 0 ); 
  output->SetSpacing( ar );
  
  memset( ar, 0, sizeof( ar ) );
  output->SetOrigin( ar );

  vtkDataArray *data = NULL;

  void *volumeDataPtr = volumeData->GetDataPtr();
  switch ( volumeData->GetType() ) 
    {
    case TYPE_BYTE: 
    {
    vtkUnsignedCharArray *byteData = vtkUnsignedCharArray::New();
    byteData->SetVoidArray( volumeDataPtr, this->m_Volume->GetNumberOfPixels(), 1 /* save */ );
    data = byteData;
    break;
    }
    case TYPE_CHAR:
    {
    vtkCharArray *charData = vtkCharArray::New();
    charData->SetVoidArray( volumeDataPtr, this->m_Volume->GetNumberOfPixels(), 1  /* save */ );
    data = charData;
    break;
    }
    case TYPE_SHORT:
    {
    vtkShortArray *shortData = vtkShortArray::New();
    shortData->SetVoidArray( volumeDataPtr, this->m_Volume->GetNumberOfPixels(), 1  /* save */ );
    data = shortData;
    break;
    }
    case TYPE_USHORT: 
    {
    vtkUnsignedShortArray *unsignedShortData = vtkUnsignedShortArray::New();
    unsignedShortData->SetVoidArray( volumeDataPtr, this->m_Volume->GetNumberOfPixels(), 1  /* save */ );
    data = unsignedShortData;
    break;
    }
    case TYPE_INT:
    {
    vtkIntArray *intData = vtkIntArray::New();
    intData->SetVoidArray( volumeDataPtr, this->m_Volume->GetNumberOfPixels(), 1  /* save */ );
    data = intData;
    break;
    }
    case TYPE_FLOAT: 
    {
    vtkFloatArray *floatData = vtkFloatArray::New();
    floatData->SetVoidArray( volumeDataPtr, this->m_Volume->GetNumberOfPixels(), 1  /* save */ );
    data = floatData;
    break;
    }
    case TYPE_DOUBLE: 
    {
    vtkDoubleArray *doubleData = vtkDoubleArray::New();
    doubleData->SetVoidArray( volumeDataPtr, this->m_Volume->GetNumberOfPixels(), 1  /* save */ );
    data = doubleData;
    break;
    }
    default:
    {
    data = NULL;
    break;
    }
    }
  
  vtkPointData* a = output->GetPointData();
  if ( data )
    {
    a->SetScalars( data );
    data->Delete();
    }
}

} // namespace cmtk
