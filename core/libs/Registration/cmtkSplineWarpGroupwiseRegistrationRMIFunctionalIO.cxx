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

#include <cmtkSplineWarpGroupwiseRegistrationRMIFunctional.h>

#include <cmtkUniformVolume.h>
#include <cmtkVolumeIO.h>
#include <cmtkClassStreamAffineXform.h>
#include <cmtkConsole.h>

#ifdef CMTK_BUILD_MPI
#  include <cmtkMPI.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ClassStream& 
operator<<
  ( ClassStream& stream, const SplineWarpGroupwiseRegistrationRMIFunctional& func )
{
  stream.Begin( "template" );
  stream.WriteIntArray( "dims", func.m_TemplateGrid->GetDims(), 3 );
  stream.WriteCoordinateArray( "delta", func.m_TemplateGrid->GetDelta(), 3 );
  stream.WriteCoordinateArray( "size", func.m_TemplateGrid->Size, 3 );
  stream.WriteCoordinateArray( "origin", func.m_TemplateGrid->m_Origin.XYZ, 3 );
  stream.End();
  
  for ( size_t idx = 0; idx < func.m_XformVector.size(); ++idx )
    {
    stream.WriteString( "target", func.m_OriginalImageVector[idx]->m_MetaInformation[CMTK_META_FS_PATH].c_str() );
    stream << func.GetXformByIndex(idx);
    }
  
  return stream;
}

ClassStream& 
operator>>
  ( ClassStream& stream, SplineWarpGroupwiseRegistrationRMIFunctional& func )
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

  UniformVolume::SmartPtr templateGrid( new UniformVolume( dims, size ) );
  templateGrid->SetOrigin( Vector3D( origin ) );

  std::vector<UniformVolume::SmartPtr> imageVector;
  std::vector<AffineXform::SmartPtr> xformVector;

  char* targetPath = stream.ReadString( "target", NULL /*default*/, false /*forward*/ );
  while ( targetPath )
    {
#ifdef CMTK_BUILD_MPI
    UniformVolume::SmartPtr image( NULL );
    if ( MPI::COMM_WORLD.Get_rank() == (imageVector.size() % MPI::COMM_WORLD.Get_size() ) )
      {
      image = UniformVolume::SmartPtr( VolumeIO::ReadOriented( targetPath, true ) );
      if ( ! image || ! image->GetData() )
	{
	StdErr << "Could not open image " << targetPath << "\n";
	exit( 1 );
	}
      }
    else
      {
      image = UniformVolume::SmartPtr( VolumeIO::ReadGridOriented( targetPath, true ) );
      }
#else
    UniformVolume::SmartPtr image( VolumeIO::ReadOriented( targetPath, true ) );
    if ( ! image || ! image->GetData() )
      {
      StdErr << "Could not open image " << targetPath << "\n";
      exit( 1 );
      }
#endif
    imageVector.push_back( image );

    AffineXform::SmartPtr xform;
    stream >> xform;
    xformVector.push_back( xform );

    free( targetPath );
    targetPath = stream.ReadString( "target", NULL /*default*/, true /*forward*/ );
    }

  func.m_InitialAffineXformsVector = xformVector;
  func.SetTargetImages( imageVector );
  func.SetTemplateGrid( templateGrid );

  return stream;
}

} // namespace cmtk
