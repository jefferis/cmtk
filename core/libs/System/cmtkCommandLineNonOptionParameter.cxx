#include <cmtkCommandLine.h>

#include <sstream>

void
cmtk::CommandLine::NonOptionParameter
::Evaluate( const size_t argc, const char* argv[], size_t& index )
{
  if ( this->Flag ) 
    *this->Flag = true;
  
  if ( index < argc ) 
    {
    *this->Var = argv[index];
    } 
  else
    {
    throw( Exception( "Argument missing", index ) );
    }
}


mxml_node_t* 
cmtk::CommandLine::NonOptionParameter
::MakeXML( mxml_node_t *const parent, const int index ) const
{
  mxml_node_t *node;
  if ( this->m_Properties & PROPS_IMAGE )
    {
    node = mxmlNewElement( parent, "image" );
    mxmlElementSetAttr( node, "type", "scalar" );
    mxmlElementSetAttr( node, "coordinateSystem", "ras" );
    }
  else if ( this->m_Properties & PROPS_FILENAME )
    node = mxmlNewElement( parent, "file" );
  else if ( this->m_Properties & PROPS_DIRNAME )
    node = mxmlNewElement( parent, "directory" );

  if ( this->m_Properties & PROPS_OUTPUT )
    mxmlNewText( mxmlNewElement( node, "channel" ), 0, "output" );
  else
    mxmlNewText( mxmlNewElement( node, "channel" ), 0, "input" );


  mxmlNewText( mxmlNewElement( node, "name" ), 0, this->m_Name );
  mxmlNewText( mxmlNewElement( node, "description" ), 0, this->m_Comment );

  std::ostringstream strm;
  strm << index;
  mxmlNewText( mxmlNewElement( node, "index" ), 0, strm.str().c_str() );
} 
