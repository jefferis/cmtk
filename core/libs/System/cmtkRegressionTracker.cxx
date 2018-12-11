/*
//
//  Copyright 2012, 2013 SRI International
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

#include "cmtkRegressionTracker.h"

#include <stdlib.h>

cmtk::RegressionTracker::RegressionTracker()
  : m_File( NULL ),
    m_WriteFlag( false )
{
  const char *env = getenv( "CMTK_RTRACKER" );
  if ( env )
    {
    this->m_File = fopen( env, "r" );
    if ( this->m_File )
      this->m_WriteFlag = false;
    else
      {
      this->m_File = fopen( env, "w" );
      this->m_WriteFlag = true;
      }
    }
}

cmtk::RegressionTracker::~RegressionTracker()
{
  if ( this->m_File )
    {
    fclose( this->m_File );
    }
}

void
cmtk::RegressionTracker::CompareChecksum( const unsigned char *const data, size_t nBytes )
{
  unsigned int checksum = 0;
  for ( size_t n = 0; n < nBytes; ++n )
    {
    checksum = ((checksum & 255) << 24) | (checksum >> 8);
    checksum ^= data[n];
    }
  
  if ( this->m_WriteFlag )
    {
    fprintf( this->m_File, "%u\n", checksum );
    }
  else
    {
    unsigned int baseline;
    if ( 1 != fscanf( this->m_File, "%20u", &baseline ) )
      {
      this->Trap();
      }
    
    if ( checksum != baseline )
      this->Trap();
    }
}
