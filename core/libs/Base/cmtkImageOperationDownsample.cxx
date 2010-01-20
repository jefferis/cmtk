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

#include <cmtkImageOperationDownsample.h>

cmtk::UniformVolume::SmartPtr
cmtk::ImageOperationDownsample
::Apply( cmtk::UniformVolume::SmartPtr& volume )
{
  const int factors[3] = { this->m_FactorX, this->m_FactorY, this->m_FactorZ };
  return cmtk::UniformVolume::SmartPtr( volume->GetDownsampled( factors ) );
}

void
cmtk::ImageOperationDownsample
::New( const char* arg )
{
  int factorsX = 1;
  int factorsY = 1;
  int factorsZ = 1;
  
  const size_t nFactors = sscanf( arg, "%d,%d,%d", &factorsX, &factorsY, &factorsZ );
  if ( nFactors == 1 )
    {
    factorsZ = factorsY = factorsX;
    }
  else
    {
    if ( nFactors != 3 )
      {
      cmtk::StdErr << "ERROR: downsampling factors must either be three integers, x,y,z, or a single integer\n";
      exit( 1 );
      }
    }
  ImageOperation::m_ImageOperationList.push_back( SmartPtr( new ImageOperationDownsample( factorsX, factorsY, factorsZ ) ) );
}

