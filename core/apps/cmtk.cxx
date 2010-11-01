/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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
#include "convert.cxx"
#include "describe.cxx"
#include "film.cxx"
#include "filter.cxx"
#include "gregxform.cxx"
#include "imagemath.cxx"
#include "levelset.cxx"
#include "mcaffine.cxx"
#include "mcwarp.cxx"
#include "model.cxx"
#include "overlap.cxx"
#include "reformatx.cxx"
#include "registration.cxx"
#include "reorient.cxx"
#include "similarity.cxx"
#include "statistics.cxx"
#include "ttest.cxx"
#include "warp.cxx"

#ifdef CMTK_HAVE_DCMTK
#include "dcm2image.cxx"
#endif

/** Set up table of command names and main function pointers */

typedef int (*testFuncPtr)( int, char*[] );
typedef struct
{
  const char* name;
  const testFuncPtr func;
} mainNameAndFunctionPointer;

const mainNameAndFunctionPointer mainFuncTable[] =
{
  { "convert",        &cmtk::apps::convert::main },
  { "describe",       &cmtk::apps::describe::main },
  { "film",           &cmtk::apps::film::main },
  { "filter",         &cmtk::apps::filter::main },
  { "gregxform",      &cmtk::apps::gregxform::main },
  { "imagemath",      &cmtk::apps::imagemath::main },
  { "levelset",       &cmtk::apps::levelset::main },
  { "mcaffine",       &cmtk::apps::mcaffine::main },
  { "mcwarp",         &cmtk::apps::mcwarp::main },
  { "model",          &cmtk::apps::model::main },
  { "overlap",        &cmtk::apps::overlap::main },
  { "reformatx",      &cmtk::apps::reformatx::main },
  { "registration",   &cmtk::apps::registration::main },
  { "reorient",       &cmtk::apps::reorient::main },
  { "similarity",     &cmtk::apps::similarity::main },
  { "statistics",     &cmtk::apps::statistics::main },
  { "ttest",          &cmtk::apps::ttest::main },
  { "warp",           &cmtk::apps::warp::main },
#ifdef CMTK_HAVE_DCMTK
  { "dcm2image",      &cmtk::apps::dcm2image::main },
#endif
  { NULL, NULL }
};

int
main( int argc, char* argv[] )
{  
  for ( const mainNameAndFunctionPointer* table = mainFuncTable; table->name; ++table )
    {
    if ( !strcmp( argv[1], table->name ) )
      return table->func( argc-1, argv+1 );
    }

  cmtk::StdErr << "ERROR: command '" << argv[1] << "' is not supported\n";
  return 1;
}
