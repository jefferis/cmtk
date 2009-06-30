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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkCommandLine_h_included_
#define __cmtkCommandLine_h_included_

#include <cmtkconfig.h>

#include <map>
#include <list>
#include <string>
#include <sstream>

#include <stdlib.h>
#include <string.h>

#include <cmtkSmartPtr.h>
#include <cmtkConsole.h>
#include <cmtkCommandLineTypeTraits.h>

#include <mxml.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Command line argument parser.
 * This class provides functionality for command line argument parsing, including automatic
 * generation of help texts and XML self-description according to the specification of the
 * Slicer3 execution model 
 * (http://www.slicer.org/slicerWiki/index.php/Slicer3:Execution_Model_Documentation)
 *
 * The class handles the following types of command line arguments:
 * 
 * 1. Switches: a switch does not have an argument itself but internally sets a variable
 *  to a pre-defined value.
 *
 * Example:
 *
 * \code
 * bool verbose;
 * cl.AddSwitch( CommandLine::Key( "verbose" ), &verbose, true, "Verbose mode );
 * \endcode
 *
 * Effect: "--verbose" sets the "verbose" variable to "true".
 *
 * 2. Options: an option has an argument, which is evaluated and stored in a variable.
 *
 * Example:
 *
 * \code
 * int iterations = 0;
 * cl.AddOption( CommandLine::Key( "iterations" ), &iterations, "Number of iterations" );
 * \endcode
 *
 * Effect: "--iterations 10" sets "iterations" to 10.
 *
 * 3. Callbacks: a callback is linked to a user-defined C function, which is called when the
 *  associated key appears on the command line.
 *
 * Example:
 *
 * \code
 * const char* callback( const char* arg )
 * {
 *    std::cerr << "callback!" << std::endl;
 * }
 * cl.AddCallback( Key( "do-something" ), &callback, "Do something using a callback function" );
 * \endcode
 */
class CommandLine
{
public:
  /// This class.
  typedef CommandLine Self;

  /// Enum for program properties.
  typedef enum
  {
    /// Program title.
    PRG_TITLE = 0,
    /// Program description.
    PRG_DESCR = 1,
    /// Program category.
    PRG_CATEG = 2,
    /// Program acknowledgement.
    PRG_ACKNL = 3,
    /// Program license (URL).
    PRG_LCNSE = 4,
    /// Program contributor.
    PRG_CNTRB = 5,
    /// Program documentation (URL).
    PRG_DOCUM = 6,
    /// Program version (default: CMTK version).
    PRG_VERSN = 7,
    /// Program syntax.
    PRG_SYNTX = 100
  } ProgramProperties;

  /// Enum for command line item properties.
  typedef enum
  {
    /// No properties.
    PROPS_NONE = 0,
    /// Item is an advanced option.
    PROPS_ADVANCED = 1,
    /// Item can appear repeatedly
    PROPS_MULTIPLE = 2,
    /// Item is not included in XML output.
    PROPS_NOXML = 4,
    /// Item is a generic file name
    PROPS_FILENAME = 8,
    /// Item is an image file name
    PROPS_IMAGE = 16,
    /// Item is an output file name
    PROPS_OUTPUT = 32
  } ItemProperties;

  /// Set program title, description, an category.
  void SetProgramInfo( const ProgramProperties key, const std::string& value )
  {
    this->m_ProgramInfo[key] = value;
  }

  /// Exception for parser error reporting.
  class Exception
  {
  public:
    /// Constructor.
    Exception( const char* message = NULL, const size_t index = 0 ) :
      Message( message ), Index( index ) 
    {}
    
    /// The error message describing the parsing condition.
    const char* Message;
    
    /// The parameter index where the condition occured.
    size_t Index;
  };

  /// Callback function.
  typedef const char* (*CallbackFunc)();

  /// Callback function with an argument.
  typedef const char* (*CallbackFuncArg)( const char* );

  /// Callback function with an integer argument.
  typedef const char* (*CallbackFuncIntArg)( const long int );

  /// Callback function with a double-precision floating point argument.
  typedef const char* (*CallbackFuncDblArg)( const double );

  /// Callback function with multiple arguments.
  typedef const char* (*CallbackFuncMultiArg)( const char**, int& argsUsed );

  /** Command line key: short and long option.
   * A key is the character (short options) or string (long options) that is associated
   * with a given action. When the key appears on the command line, the action is executed.
   *
   * The following are examples of Key initializers:
   * \code
   * Key( 'v', "verbose" );
   * Key( "number-of-iterations" );
   * Key( 'q' );
   * \endcode
   */
  class Key
  {
  public:
    /** Short option constructor.
     *\deprecated This may be removed in the future because short-only options are not nice when we use
     * the XML self description for Slicer integration.
     */
    Key( const char keyChar ) : m_KeyChar( keyChar ) 
    {
      StdErr << "WARNING: short command line option '" << keyChar << "' should also have a long name.\n";
    }

    /// Long option constructor.
    Key( const std::string& keyString ) : m_KeyChar( 0 ), m_KeyString( keyString ) {}

    /// Short and long option constructor.
    Key( const char keyChar, const std::string& keyString ) : m_KeyChar( keyChar ), m_KeyString( keyString ) {}

    /// Copy constructor.
    Key( const Key& other ) : m_KeyChar( other.m_KeyChar ), m_KeyString( other.m_KeyString ) {}

    /// Short option key.
    char m_KeyChar;

    /// Long option key.
    std::string m_KeyString;
  };

  /// Virtual base class for command line item.
  class Item
  {
  public:
    /// Smart pointer to this class.
    typedef SmartPointer<Item> SmartPtr;

    /// Constructor.
    Item() : m_Properties( PROPS_NONE ) {};

    /// Virtual destructor.
    virtual ~Item() {}

    /// Set item properties.
    virtual void SetProperties( const ItemProperties properties )
    {
      this->m_Properties = properties;
    }
    
    /// Virtual function: evaluate switch or option.
    virtual void Evaluate( const size_t argc, const char* argv[], size_t& index ) = 0;

    /// Virtual function that returns an XML tree describing this option.
    virtual mxml_node_t* MakeXML( mxml_node_t *const parent //!< Parent in the XML tree for the new node.
      ) const = 0;    

  protected:
    /// Item properties.
    ItemProperties m_Properties;

    /// Safely convertstring argument to long integer.
    static long int ConvertStrToLong( const char* str );

    /// Safely convert string argument to double.
    static double ConvertStrToDouble( const char* str );

    /// Convert function.
    template<class T> T Convert( const char* str );

  private:
    /// Allow command line class full access.
    friend class CommandLine;
  };

private:
  /// Command line switch.
  template<class T> 
  class Switch : 
    /// Inherit from generic command line item.
    public Item
  {
  public:
    /// Constructor.
    Switch( T *const flag, const T value ) : Flag( flag ), Value( value ) {}

    /// Evaluate and set associated flag.
    virtual void Evaluate( const size_t, const char*[], size_t& )
    { 
      *Flag  = Value;
    }

    /// Virtual function that returns an XML tree describing this option.
    virtual mxml_node_t* MakeXML(  mxml_node_t *const parent ) const 
    {
      mxml_node_t *node = mxmlNewElement( parent, "boolean" );
      return node;
    }

  private:
    /// Pointer to flag handled by this switch.
    T* Flag;

    /// Value to set flag to when switch is encountered.
    const T Value;
  };

  /// Command line option with argument.
  template<class T>
  class Option : 
    /// Inherit from generic command line item.
    public Item
  {
  public:
    /// Constructor.
    Option( T *const var, bool *const flag ) : Var( var ), Flag( flag ) {}
    
    /// Evaluate and set associated option.
    virtual void Evaluate( const size_t argc, const char* argv[], size_t& index )
    {
      if ( Flag ) *Flag = true;
      if ( index+1 < argc ) 
	{
	*Var = this->Convert<T>( argv[index+1] );
	++index;
	} 
      else
	{
	throw( Exception( "Option needs an argument.", index ) );
	}
    }

    /// Virtual function that returns an XML tree describing this option.
    virtual mxml_node_t* MakeXML(  mxml_node_t *const parent ) const 
    {
      mxml_node_t *node = mxmlNewElement( parent, CommandLineTypeTraits<T>::GetName() );
      if ( !Flag ) // if there is no flag monitoring this option, then there must be a valid default value
	{
	mxml_node_t *dflt = mxmlNewElement( node, "default" );
	mxmlNewText( dflt, 0, CommandLineTypeTraits<T>::ValueToString( *Var ).c_str() );
	}
      return node;
    }

  private:
    /// Pointer to associated variable.
    T* Var;

    /// Pointer to (optional) flag that will be set when this option occurs.
    bool* Flag;
  };

  /// Command line option with argument.
  template<class T>
  class Repeat : 
    /// Inherit from generic command line item.
    public Item
  {
  public:
    /// Constructor.
    Repeat( std::list<T>& list ) : m_pList( &list ) 
    {
      if ( ! m_pList )
	throw( Exception( "List pointer must not be NULL" ) );
    }
    
    /// Evaluate and set associated option.
    virtual void Evaluate( const size_t argc, const char* argv[], size_t& index )
    {
      if ( index+1 < argc ) 
	{
	m_pList->push_back( this->Convert<T>( argv[index+1] ) );
	++index;
	} 
      else
	{
	throw( Exception( "Option needs an argument.", index ) );
	}
    }

    /// Virtual function that returns an XML tree describing this option.
    virtual mxml_node_t* MakeXML(  mxml_node_t *const parent ) const 
    {
      mxml_node_t *node = mxmlNewElement( parent, CommandLineTypeTraits<T>::GetName() );
      mxmlNewText( mxmlNewElement( node, "name" ), 0, CommandLineTypeTraits<T>::GetName() );
      return node;
    }

  private:
    /// Pointer to associated variable.
    std::list<T>* m_pList;
  };

  /// Command line callback.
  class Callback : 
    /// Inherit from generic command line item.
    public Item
  {
  public:
    /// Constructor for simple callback.
    Callback( CallbackFunc func ) : m_Func( func ), m_FuncArg( NULL ), m_FuncIntArg( NULL ), m_FuncDblArg( NULL ), m_FuncMultiArg( NULL ) {}
    
    /// Constructor for callback with argument.
    Callback( CallbackFuncArg funcArg ) : m_Func( NULL ), m_FuncArg( funcArg ), m_FuncIntArg( NULL ), m_FuncDblArg( NULL ), m_FuncMultiArg( NULL ) {}
    
    /// Constructor for callback with integer argument.
    Callback( CallbackFuncIntArg funcIntArg ) : m_Func( NULL ), m_FuncArg( NULL ), m_FuncIntArg( funcIntArg ), m_FuncDblArg( NULL ), m_FuncMultiArg( NULL ) {}
    
    /// Constructor for callback with integer argument.
    Callback( CallbackFuncDblArg funcDblArg ) : m_Func( NULL ), m_FuncArg( NULL ), m_FuncIntArg( NULL ), m_FuncDblArg( funcDblArg ), m_FuncMultiArg( NULL ) {}
    
    /// Constructor for callback with multiple arguments.
    Callback( CallbackFuncMultiArg funcMultiArg ) : m_Func( NULL ), m_FuncArg( NULL ), m_FuncIntArg( NULL ), m_FuncDblArg( NULL ), m_FuncMultiArg(  funcMultiArg ) {}
    
    /// Evaluate and set associated option.
    virtual void Evaluate( const size_t argc, const char* argv[], size_t& index );

    /// Virtual function that returns an XML tree describing this option.
    virtual mxml_node_t* MakeXML(  mxml_node_t *const parent ) const;

  private:
    /// Pointer to simple callback.
    CallbackFunc m_Func;

    /// Pointer to callback with argument.
    CallbackFuncArg m_FuncArg;

    /// Pointer to callback with integer argument.
    CallbackFuncIntArg m_FuncIntArg;

    /// Pointer to callback with floating point argument.
    CallbackFuncDblArg m_FuncDblArg;

    /// Pointer to callback with multiple arguments.
    CallbackFuncMultiArg m_FuncMultiArg;
  };

public:
  /// Constructor.
  CommandLine( int argc, char* argv[], const int properties = PROPS_NONE ) 
  {
    this->SetDefaultInfo();    
    ArgC = argc;
    ArgV = const_cast<const char**>( argv );
    this->m_Properties = properties;

    this->BeginGroup( "MAIN", "Main Options" );
  }

  /// Constructor.
  CommandLine( const int argc, const char* argv[], const int properties = PROPS_NONE ) 
  {
    this->SetDefaultInfo();
    ArgC = argc;
    ArgV = argv;
    this->m_Properties = properties;

    this->BeginGroup( "MAIN", "Main Options" );
  }

  /// Destructor.
  ~CommandLine() {};

  /// Add switch.
  template<class T> 
  Item::SmartPtr& 
  AddSwitch( const Key& key, T *const var, const T value, const char* comment ) 
  {
    return this->AddKeyAction( KeyToAction::SmartPtr( new KeyToAction( key, Item::SmartPtr( new Switch<T>( var, value ) ), comment ) ) );
  }
  
  /// Add switch for const char* variables.
  Item::SmartPtr&
  AddSwitch( const Key& key, const char** var, const char* value, const char* comment ) 
  {
    return this->AddKeyAction( KeyToAction::SmartPtr( new KeyToAction( key, Item::SmartPtr( new Switch<const char*>( var, value ) ), comment ) ) );
  }

  /// Add option.
  template<class T> 
  Item::SmartPtr&
  AddOption( const Key& key, T *const var, const char* comment, bool *const flag = NULL ) 
  {
    return this->AddKeyAction( KeyToAction::SmartPtr( new KeyToAction( key, Item::SmartPtr( new Option<T>( var, flag ) ), comment ) ) );
  }
  
  /// Add repeat option (put arguments into a FIFO list).
  template<class T> 
  Item::SmartPtr&
  AddRepeat( const Key& key, std::list<T>& list, const char* comment ) 
  {
    return this->AddKeyAction( KeyToAction::SmartPtr( new KeyToAction( key, Item::SmartPtr( new Repeat<T>( list ) ), comment ) ) );
  }
  
  /// Add callback.
  Item::SmartPtr&
  AddCallback( const Key& key, CallbackFunc func, const char* comment ) 
  {
    return this->AddKeyAction( KeyToAction::SmartPtr( new KeyToAction( key, Item::SmartPtr( new Callback( func ) ), comment ) ) );
  }

  /// Add callback with a single argument.
  template<class TArg>
  Item::SmartPtr&
  AddCallback( const Key& key, const char* (*funcArg)( const TArg ), const char* comment ) 
  {
    return this->AddKeyAction( KeyToAction::SmartPtr( new KeyToAction( key, Item::SmartPtr( new Callback( funcArg ) ), comment ) ) );
  }
  
  /// Add callback with multiple arguments.
  Item::SmartPtr&
  AddCallback( const Key& key, CallbackFuncMultiArg funcMultiArg, const char* comment ) 
  {
    return this->AddKeyAction( KeyToAction::SmartPtr( new KeyToAction( key, Item::SmartPtr( new Callback( funcMultiArg ) ), comment ) ) );
  }
  
  /// Begin parameter group.
  void BeginGroup( const char* name, const char* description );

  /// End last parameter group.
  void EndGroup();

  /// Parse command line.
  bool Parse();

  /// Print help text.
  void PrintHelp() const;

  /// Get next parameter index.
  size_t GetIndex() const { return Index; }

  /** Get next command line argument.
   * An exception is generated if no further arguments are available.
   */
  const char* GetNext() 
  {
    if ( Index >= ArgC ) 
      throw( Exception( "Out of arguments.", Index ) ); 
    return ArgV[Index++];
  }
  
  /** Get next command line argument.
   * A null pointer is returned if no further arguments are available.
   */
  const char* GetNextOptional() 
  {
    if ( Index >= ArgC ) return NULL;
    return ArgV[Index++];
  }

  /** Write XML self description according to Slice3 execution model.
   *@see 
   */
  void WriteXML() const;
  
private:
  /// Total number of arguments.
  size_t ArgC;

  /// Array of argument pointers.
  const char** ArgV;

  /// Global properties of the command line.
  int m_Properties;

  /// Index of current argument.
  size_t Index;

  /// Local class that connects command line options with their evaluators.
  class KeyToAction
  {
  public:
    /// Smart pointer to this class.
    typedef SmartPointer<KeyToAction> SmartPtr;

    /// Constructor.
    KeyToAction( const Key& key, //!< Key: long and/or short command line option for this action.
		 Item::SmartPtr action, //!< The actual action (option, switch, callback, etc.)
		 const char* comment ) : //!< Command line help comment for this action.
      m_Key( key.m_KeyChar ),
      m_KeyString( key.m_KeyString ),
      m_Action( action ),
      m_Comment( comment )
    {}
    
    /// Returns an XML tree describing this key and action.
    mxml_node_t* MakeXML( mxml_node_t *const parent //!< Parent in the XML tree for the new node.
      ) const;

    /// Print help for this item.
    void PrintHelp( const size_t globalIndent = 0 ) const;
    
  private:
    /// Short option associated with this action.
    char m_Key;

    /// Long option associated with this action.
    std::string m_KeyString;
    
    /// Action.
    Item::SmartPtr m_Action;

    /// Comment (description).
    const char* m_Comment;

    /// Give command line class full access.
    friend class CommandLine;
  };

  /// Type for key/action lists.
  typedef std::list<KeyToAction::SmartPtr> KeyActionListType;

  /// Reference to current key/action list (current group).
  KeyActionListType* m_KeyActionList;
				    
  /// Reference to complete key/action list (all groups combined), which is used for key lookup.
  KeyActionListType m_KeyActionListComplete;
				    
  /// Add an item to applicable key/action lists.
  Item::SmartPtr& AddKeyAction( KeyToAction::SmartPtr keyAction )
  {
    this->m_KeyActionList->push_back( keyAction );
    this->m_KeyActionListComplete.push_back( keyAction );
    return this->m_KeyActionListComplete.back()->m_Action;
  }

  /// Type for action groups, which map a group name to the group's key-action list.
  class KeyActionGroupType
  {
  public:
    /// Smart pointer to this class.
    typedef SmartPointer<KeyActionGroupType> SmartPtr;
    
    /// Constructor: set name and description.
    KeyActionGroupType( const std::string& name, const std::string& description ) : m_Name( name ), m_Description( description ) {};

    /// Group name.
    std::string m_Name;

    /// Group description.
    std::string m_Description;

    /// Key-action list for this group.
    KeyActionListType m_KeyActionList;
  };

  /// Type for group list.
  typedef std::list<KeyActionGroupType::SmartPtr> KeyActionGroupListType;
  
  /// List of command line keys with associated actions.
  KeyActionGroupListType m_KeyActionGroupList;

  /// Map type for program meta information.
  typedef std::map<ProgramProperties,std::string> ProgramPropertiesMapType;

  /// Program meta information.
  ProgramPropertiesMapType m_ProgramInfo;

  /// Add program info item to XML tree.
  mxml_node_t* AddProgramInfoXML( mxml_node_t *const parent, //!< Parent node for new entry in XML tree.
				  const ProgramProperties key, //!< Key code for program property.
				  const char* name ) const; //!< Name of XML tag for this property.

  /// Set default values for meta information.
  void SetDefaultInfo();
  
  /// Make options class friend.
  template<class T> friend class Option;  

  /// Make switch class friend.
  template<class T> friend class Switch;  

  /// Make callback class friend.
  friend class Callback;  
};

/// Output of command line exception.
Console& operator<<( Console& console, CommandLine::Exception e );

//@}

} // namespace cmtk

#include <cmtkCommandLineConvert.txx>

#endif // #ifndef __cmtkCommandLine_h_included_
