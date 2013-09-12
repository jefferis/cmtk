/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2009, 2013 SRI International
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

#include <Registration/cmtkProtocolCallback.h>

#include <Base/cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ProtocolCallback::ProtocolCallback
( const std::string& filename, const bool debug ) 
{ 
  if ( !filename.empty() ) 
    {
    if ( (fp = fopen( filename.c_str(), "w" )) ) 
      {
      fputs( "4\n1 3 3 3\n", fp );
      fflush( fp );
      }
    }
  else
    fp = NULL; 
  
  Debug = debug;
}

ProtocolCallback::~ProtocolCallback () 
{
  if (fp) fclose(fp); 
}

CallbackResult
ProtocolCallback::ExecuteWithData
( const CoordinateVector& v, const double metric ) 
{
  size_t dim = std::min<unsigned int>( 20, v.Dim );
  if (fp) 
    {
    fprintf( fp, "%f", metric );    
    for ( size_t i = 0; i < dim; ++i )
      fprintf( fp, " %f", (float) v[i] );
    
    if ( v.Dim > 20 ) fputs( " ...", fp );
    fputs( "\n", fp );
    fflush( fp );
    }
  
  if ( Debug ) 
    {
    fprintf( stderr, "%f", metric );
    for ( size_t i = 0; i < dim; ++i )
      fprintf( stderr, " %f", (float) v[i] );
    fputs( "\n", stderr );
    }
  
  return this->Superclass::ExecuteWithData( v, metric );
}

void
ProtocolCallback::Comment ( const std::string& comment )
{
  if ( fp ) 
    {
    if ( !comment.empty() ) 
      {
      fprintf( fp, "# %s\n", comment.c_str() );
      fflush( fp );
      } 
    else
      {
      fputs( "#\n", fp );
      fflush( fp );
      }
    }
  
  if ( Debug )
    {
    if ( !comment.empty() )
      fprintf( stderr, "# %s\n", comment.c_str() );
    else
      fputs( "#\n", stderr );
    }
}

} // namespace cmtk
