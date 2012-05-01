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
    { "phantomType",  { NULL, NULL, NULL, "\n" } },
    { "snr",          { NULL, NULL, NULL, "\n" } },
    { "cnr",          { NULL, NULL, NULL, "\n" } },
    { "landmarkList", { NULL, "\n", NULL, "\n" } },
    { "landmark",     { "\t", "\n", "\t", "\n" } },
    { "name",         { "\t\t", NULL, NULL, "\n" } },
    { "detected",     { "\t\t", NULL, NULL, "\n" } },
    { "expected",     { "\t\t", NULL, NULL, "\n" } },
    { "isPrecise",    { "\t\t", NULL, NULL, "\n" } },
    { "fitResidual",  { "\t\t", NULL, NULL, "\n" } },
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
  mxmlNewText( mxmlNewElement( x_phantom, "phantomType" ), 0, "MagphanEMR051" );
  mxmlNewReal( mxmlNewElement( x_phantom, "snr" ), phantom.m_EstimatedSNR );    
  mxmlNewReal( mxmlNewElement( x_phantom, "cnr" ), phantom.m_EstimatedCNR );    
    
  mxml_node_t *x_lmpairs = mxmlNewElement( x_phantom, "landmarkList" );
  mxmlElementSetAttr( x_lmpairs, "coordinates", "physical" );
  mxmlElementSetAttr( x_lmpairs, "space", "RAS" );

  const std::list<LandmarkPair>& lmPairs = phantom.LandmarkPairsList();
  for ( std::list<LandmarkPair>::const_iterator it = lmPairs.begin(); it != lmPairs.end(); ++it )
    {
    mxml_node_t *x_lm = mxmlNewElement( x_lmpairs, "landmark");
    
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
    
    mxmlNewText( mxmlNewElement( x_lm, "isPrecise" ), 0, it->m_Precise ? "yes" : "no" );
    mxmlNewReal( mxmlNewElement( x_lm, "fitResidual" ), it->m_Residual );    
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

cmtk::DetectedPhantomMagphanEMR051::SmartPtr
cmtk::PhantomIO::Read( const std::string& fpath )
{
  FILE *file = fopen( fpath.c_str(), "r" );
  if ( !file )
    {
    cmtk::StdErr << "ERROR: could not open file " << fpath << " for reading\n";
    return DetectedPhantomMagphanEMR051::SmartPtr( NULL );
    }

  mxml_node_t *x_root = mxmlLoadFile( NULL, file, MXML_TEXT_CALLBACK );
  fclose( file );

  mxml_node_t *x_landmarks = mxmlFindElement( x_root, x_root, "landmarkList", NULL, NULL, MXML_DESCEND );
  if ( ! x_landmarks )
    {
    cmtk::StdErr << "ERROR: could not file 'landmarks' XML element in file " << fpath << "\n";
    mxmlDelete( x_root );
    return DetectedPhantomMagphanEMR051::SmartPtr( NULL );
    }

  AffineXform xform;
  DetectedPhantomMagphanEMR051::SmartPtr result( new DetectedPhantomMagphanEMR051( 0, 0, xform ) );
  
  for ( mxml_node_t* x_fiducial = mxmlFindElement( x_landmarks, x_root, "landmark", NULL, NULL, MXML_DESCEND ); x_fiducial != NULL; x_fiducial = mxmlFindElement( x_fiducial, x_root, "landmark", NULL, NULL, MXML_DESCEND ) )
    {
    mxml_node_t* x_name = mxmlFindElement( x_fiducial, x_root, "name", NULL, NULL, MXML_DESCEND );
    if ( ! x_name || ! x_name->child ) continue;
    const std::string name = x_name->child->value.text.string;

    mxml_node_t* x_expected = mxmlFindElement( x_fiducial, x_root, "expected", NULL, NULL, MXML_DESCEND );
    if ( ! x_expected || ! x_expected->child ) continue;

    Landmark::SpaceVectorType expected;
    mxml_node_t* x_expected_it = x_expected->child;
    for ( size_t i = 0; i < 3; ++i, x_expected_it = x_expected_it->next )
      expected[i] = atof( x_expected_it->value.text.string );

    mxml_node_t* x_detected = mxmlFindElement( x_fiducial, x_root, "detected", NULL, NULL, MXML_DESCEND );
    if ( ! x_detected || ! x_detected->child ) continue;

    Landmark::SpaceVectorType detected;
    mxml_node_t* x_detected_it = x_detected->child;
    for ( size_t i = 0; i < 3; ++i, x_detected_it = x_detected_it->next )
      detected[i] = atof( x_detected_it->value.text.string );
    
    mxml_node_t* x_precise = mxmlFindElement( x_fiducial, x_root, "isPrecise", NULL, NULL, MXML_DESCEND );
    if ( ! x_precise || ! x_precise->child ) continue;
    const bool precise = !strcmp( x_precise->child->value.text.string, "yes" );

    mxml_node_t* x_residual = mxmlFindElement( x_fiducial, x_root, "fitResidual", NULL, NULL, MXML_DESCEND );
    if ( ! x_residual || ! x_residual->child ) continue;
    const Types::Coordinate residual = static_cast<Types::Coordinate>( atof( x_residual->child->value.text.string ) );

    result->AddLandmarkPair( name, expected, detected, residual, precise );
    }
  
  mxmlDelete( x_root );
  
  return result;
}
