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

#include <Registration/cmtkGroupwiseRegistrationFunctionalXformTemplate.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkClassStreamAffineXform.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ClassStreamOutput& 
operator<<
  ( ClassStreamOutput& stream, const GroupwiseRegistrationFunctionalXformTemplate<AffineXform>& func )
{
  stream.Begin( "template" );
  stream.WriteIntArray( "dims", func.m_TemplateGrid->GetDims().begin(), 3 );
  stream.WriteCoordinateArray( "delta", func.m_TemplateGrid->Deltas().begin(), 3 );
  stream.WriteCoordinateArray( "size", func.m_TemplateGrid->m_Size.begin(), 3 );
  stream.WriteCoordinateArray( "origin", func.m_TemplateGrid->m_Offset.begin(), 3 );
  stream.End();
  
  for ( size_t idx = 0; idx < func.m_XformVector.size(); ++idx )
    {
    stream.WriteString( "target", func.m_OriginalImageVector[idx]->GetMetaInfo( META_FS_PATH ).c_str() );
    stream << (*func.GetXformByIndex( idx ));
    }
  
  return stream;
}

ClassStreamInput& 
operator>>
  ( ClassStreamInput& stream, GroupwiseRegistrationFunctionalXformTemplate<AffineXform>& func )
{
  if ( ! stream.Seek( "template" ) )
    {
    StdErr << "ERROR: no 'template' section in input archive\n";
    return stream;
    }

  int dims[3];
  stream.ReadIntArray( "dims", dims, 3 );
  Types::Coordinate size[3];
  stream.ReadCoordinateArray( "size", size, 3 );
  Types::Coordinate origin[3];
  stream.ReadCoordinateArray( "origin", origin, 3 );
  stream.End();

  UniformVolume::SmartPtr templateGrid( new UniformVolume( UniformVolume::IndexType::FromPointer( dims ), UniformVolume::CoordinateVectorType::FromPointer( size ) ) );
  templateGrid->SetOffset( FixedVector<3,Types::Coordinate>::FromPointer( origin ) );

  std::vector<UniformVolume::SmartPtr> imageVector;
  std::vector<AffineXform::SmartPtr> xformVector;

  char* targetPath = stream.ReadString( "target", NULL /*default*/, false /*forward*/ );
  while ( targetPath )
    {
    UniformVolume::SmartPtr image( VolumeIO::ReadOriented( targetPath ) );
    if ( ! image || ! image->GetData() )
      {
      StdErr << "ERROR: Could not open image " << targetPath << "\n";
      exit( 1 );
      }
    imageVector.push_back( image );

    AffineXform::SmartPtr xform;
    stream >> xform;
    xformVector.push_back( xform );

    free( targetPath );
    targetPath = stream.ReadString( "target", NULL /*default*/, true /*forward*/ );
    }

  func.SetTargetImages( imageVector );
  func.SetTemplateGrid( templateGrid );
  func.SetXforms( xformVector );

  return stream;
}

} // namespace cmtk
