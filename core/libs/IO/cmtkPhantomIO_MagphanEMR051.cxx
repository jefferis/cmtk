/*
//
//  Copyright 2012 SRI International
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

#include "cmtkPhantomIO.h"

#include <stdio.h>

const char *
cmtk::PhantomIO::WhitespaceWriteMiniXML( mxml_node_t* node, int where)
{
  const char* name = node->value.element.name;
  
  typedef struct _wsLookupType
  {
    /// XML element name.
    const char* name;
    /// Table of whitespace sequences.
    const char* ws[4];
  } wsLookupType;

  static const wsLookupType wsLookup[] = 
  {
    { "type",   { NULL, NULL, NULL, "\n" } },
    { "snr",       { NULL, NULL, NULL, "\n" } },
    { "cnr",       { NULL, NULL, NULL, "\n" } },
    { "landmarks", { NULL, "\n", NULL, "\n" } },
    { "fiducial",  { "\t", "\n", "\t", "\n" } },
    { "name",      { "\t\t", NULL, NULL, "\n" } },
    { "detected",  { "\t\t", NULL, NULL, "\n" } },
    { "expected",  { "\t\t", NULL, NULL, "\n" } },
    { "precise",   { "\t\t", NULL, NULL, "\n" } },
    { "residual",  { "\t\t", NULL, NULL, "\n" } },
    { NULL, {NULL, NULL, NULL, NULL} }
  };
  
  if ( (where >= 0) && (where < 4) )
    {
    for ( size_t idx = 0; wsLookup[idx].name; ++idx )
      {
      if ( ! strcmp( name, wsLookup[idx].name ) )
	return wsLookup[idx].ws[where];
      }
    }

  switch ( where )
    {
    case MXML_WS_BEFORE_OPEN:
      return NULL;
    case MXML_WS_AFTER_OPEN:
      return "\n";
    case MXML_WS_BEFORE_CLOSE:
      return NULL;
    case MXML_WS_AFTER_CLOSE:
      return "\n";
    }

  return NULL;
}

void
cmtk::PhantomIO::Write( const DetectedPhantomMagphanEMR051& phantom, const std::string& fpath )
{
  mxmlSetWrapMargin( 120 ); // make enough room for indented landmark locations
  mxml_node_t *x_root = mxmlNewElement( NULL, "?xml version=\"1.0\" encoding=\"utf-8\"?" );
  
  mxml_node_t *x_phantom = mxmlNewElement( x_root, "phantom" );
  mxmlNewText( mxmlNewElement( x_phantom, "type" ), 0, "MagphanEMR051" );
  mxmlNewReal( mxmlNewElement( x_phantom, "snr" ), phantom.m_EstimatedSNR );    
  mxmlNewReal( mxmlNewElement( x_phantom, "cnr" ), phantom.m_EstimatedCNR );    
    
  mxml_node_t *x_lmpairs = mxmlNewElement( x_phantom, "landmarks" );
  mxmlElementSetAttr( x_lmpairs, "coordinates", "Physical" );
  mxmlElementSetAttr( x_lmpairs, "space", "RAS" );

  const std::list<LandmarkPair>& lmPairs = phantom.LandmarkPairsList();
  for ( std::list<LandmarkPair>::const_iterator it = lmPairs.begin(); it != lmPairs.end(); ++it )
    {
    mxml_node_t *x_lm = mxmlNewElement( x_lmpairs, "fiducial");
    
    mxmlNewText( mxmlNewElement( x_lm, "name" ), 0, it->m_Name.c_str() ); 
    mxml_node_t *x_expected = mxmlNewElement( x_lm, "expected");
    for ( size_t idx = 0; idx < 3; ++idx )
      {
      mxmlNewReal( x_expected, it->m_Location[idx] );
      }
    
    mxml_node_t *x_detected= mxmlNewElement( x_lm, "detected");
    for ( size_t idx = 0; idx < 3; ++idx )
      {
      mxmlNewReal( x_detected, it->m_TargetLocation[idx] );
      }
    
    mxmlNewText( mxmlNewElement( x_lm, "precise" ), 0, it->m_Precise ? "yes" : "no" );
    mxmlNewReal( mxmlNewElement( x_lm, "residual" ), it->m_Residual );    
    }

  FILE *file = fopen( fpath.c_str(), "w" );
  if ( file )
    {
    mxmlSaveFile( x_root, file, Self::WhitespaceWriteMiniXML );
    fputs( "\n", file ); // end last line
    fclose( file );
    }
  else
    {
    cmtk::StdErr << "ERROR: could not open file " << fpath << " for writing\n";
    }
  
  mxmlDelete( x_root );
  
}
