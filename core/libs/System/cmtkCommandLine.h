/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <System/cmtkCannotBeCopied.h>

#include <map>
#include <list>
#include <string>
#include <sstream>
#include <vector>

#include <stdlib.h>
#include <string.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkConsole.h>
#include <System/cmtkCommandLineTypeTraits.h>
#include <System/cmtkExitException.h>

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
 *\section secSwitch Switches
 *
 *  A command line switch does not have an argument itself but internally sets a variable
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
 *\section secOption Options
 *
 * A command line option has an argument, which is evaluated and stored in a variable.
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
 *\section secCallback Callbacks
 *
 *  A callback is linked to a user-defined C function, which is called when the
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
 *
 *\section secEnumeration Enumerations
 *
 * An enumeration is a group of options that modify the same variable by setting it to
 *  different values.
 *
 * Examples:
 *
 * \code
 *  int InterleaveAxis = 1;
 *  cmtk::CommandLine::EnumGroup<int>::SmartPtr interleaveGroup = cl.AddEnum( "interleave-axis", &InterleaveAxis, "Define interleave axis." );
 *  interleaveGroup->AddSwitch( Key( "guess-from-input" ), -1, "Guess from input image" );
 *  interleaveGroup->AddSwitch( Key( 'a', "axial" ), 2, "Interleaved axial images" );
 *  interleaveGroup->AddSwitch( Key( 'c', "coronal" ), 1, "Interleaved coronal images" );
 *  interleaveGroup->AddSwitch( Key( 's', "sagittal" ), 0, "Interleaved sagittal images" );
 * \endcode
 *
 * \code
 *  std::string channel( "spgr" );
 *  cmtk::CommandLine::EnumGroup<std::string>::SmartPtr channelGroup = cl.AddEnum( "RegistrationChannel", &channel, "MR channel." );
 *  channelGroup->AddSwitch( Key( "spgr" ), "spgr", "SPGR (T1-weighted) structural channel" );
 *  channelGroup->AddSwitch( Key( "early-fse" ), "erly", "Early-echo (PD-weighted) fast spin echo channel" );
 *  channelGroup->AddSwitch( Key( "late-fse" ), "late", "Late-echo (T2-weighted) fast spin echo channel" );
 * \endcode
 */
class CommandLine :
  /// Make class uncopyable via inheritance.
  private CannotBeCopied
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
    /// No properties, but supports XML (NOXML is off)..
    PROPS_XML = 0,
    /// Item is an advanced option.
    PROPS_ADVANCED = 1,
    /// Item can appear repeatedly
    PROPS_MULTIPLE = 2,
    /// Item is not included in XML output.
    PROPS_NOXML = 4,
    /// Item is a generic file name
    PROPS_DIRNAME = 8,
    /// Item is a generic directory name
    PROPS_FILENAME = 16,
    /// Item is an image file name
    PROPS_IMAGE = 32,
    /// When used with PROPS_IMAGE, this means the expected image is a label map.
    PROPS_LABELS = 64,
    /// Item is a transformation file name
    PROPS_XFORM = 128,
    /// This parameter refers to an output, not an input.
    PROPS_OUTPUT = 256,
    /// This parameter (non-option argument) is optional.
    PROPS_OPTIONAL = 512
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
    
    /// Constructor.
    Exception( const std::string& message, const size_t index = 0 ) :
      Message( message ), Index( index ) 
    {}
    
    /// The error message describing the parsing condition.
    std::string Message;
    
    /// The parameter index where the condition occured.
    size_t Index;
  };

  /// Callback function.
  typedef void (*CallbackFunc)();

  /// Callback function with an argument.
  typedef void (*CallbackFuncArg)( const char* );

  /// Callback function with an integer argument.
  typedef void (*CallbackFuncIntArg)( const long int );

  /// Callback function with a double-precision floating point argument.
  typedef void (*CallbackFuncDblArg)( const double );

  /// Callback function with multiple arguments.
  typedef void (*CallbackFuncMultiArg)( const char**, int& argsUsed );

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
    /// Long option constructor.
    explicit Key( const std::string& keyString ) : m_KeyChar( 0 ), m_KeyString( keyString ) {}

    /// Short and long option constructor.
    explicit Key( const char keyChar, const std::string& keyString ) : m_KeyChar( keyChar ), m_KeyString( keyString ) {}

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
    Item() : m_Properties( PROPS_XML ) {};

    /// Virtual destructor.
    virtual ~Item() {}

    /// Set item properties.
    virtual Item* SetProperties( const long int properties )
    {
      this->m_Properties = properties;
      return this;
    }

    /// Get item properties.
    virtual long int GetProperties() const
    {
      return this->m_Properties;
    }

    /// Set an attribute.
    virtual Item* SetAttribute( const std::string& key, const std::string& value )
    {
      this->m_Attributes[key] = value;
      return this;
    }
    
    /// Virtual function: evaluate switch or option.
    virtual void Evaluate( const size_t argc, const char* argv[], size_t& index ) = 0;

    /// Virtual function that returns an XML tree describing this option.
    virtual mxml_node_t* MakeXML( mxml_node_t *const parent /*!< Parent in the XML tree for the new node.*/ ) const = 0;

    /// Virtual function returns a string that describes the parameter type associated with this option (derived classes only)
    virtual std::string GetParamTypeString() const { return ""; }

    /// Format additional help information.
    virtual std::ostringstream& PrintHelp( std::ostringstream& fmt /*!< Stream that the additional help information is formatted into*/ ) const
    {
      // by default, simply return stream unchanged.
      return fmt;
    }

    /// Format additional help information.
    virtual void PrintWiki() const
    {
      // by default, do nothing
    }

    /// Return true if and only if this item is the default for the associated action or variable.
    virtual bool IsDefault() const
    {
      return false;
    }

  protected:
    /// Item properties.
    long int m_Properties;

    /// Item attributes. These are free-form string key/value pairs.
    std::map<std::string, std::string> m_Attributes;

    /// Safely convertstring argument to long integer.
    static long int ConvertStrToLong( const char* str );

    /// Safely convert string argument to double.
    static double ConvertStrToDouble( const char* str );

    /// Convert function.
    template<class T> T Convert( const char* str );

    /// Helper class to avoid function template
    template<class T>
    class Helper
    {
    public:
      /// Make a basic XML node for a generic item based on type, properties, attribute, etc.
      static mxml_node_t* MakeXML( const Item* item, mxml_node_t *const parent );

      /// Return a textual description of the parameter associated with this option
      static std::string GetParamTypeString( const Item* item );
    };

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
    Switch( T *const field, const T value ) : 
      m_Field( field ), 
      m_Value( value ) 
    {
    }

    /// Virtual destructor.
    virtual ~Switch() {}

    /// Evaluate and set associated flag.
    virtual void Evaluate( const size_t, const char*[], size_t& )
    { 
      *this->m_Field  = this->m_Value;
    }

    /// Virtual function that returns an XML tree describing this option.
    virtual mxml_node_t* MakeXML(  mxml_node_t *const parent ) const 
    {
      if ( ! (this->m_Properties & PROPS_NOXML) )
	{
	return mxmlNewElement( parent, "boolean" );
	}
      return NULL;
    }

    /// Format additional help information (e.g., default values).
    virtual std::ostringstream& PrintHelp( std::ostringstream& fmt /*!< Stream that the additional help information is formatted into*/ ) const
    {
      if ( this->IsDefault() )
	fmt << "\n[This is the default]";
      return fmt;
    }
    
    /// Format additional help information (e.g., default values).
    virtual void PrintWiki() const
    {
      if ( this->IsDefault() )
	StdOut << " '''[This is the default]'''";
    }
    
    /// Return true if and only if this item is the default for the associated action or variable.
    virtual bool IsDefault() const
    {
      return ( *(this->m_Field) == this->m_Value );
    }

  private:
    /// Pointer to field handled by this switch.
    T* m_Field;

    /// Value to set field to when switch is encountered.
    const T m_Value;
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
    
    /// Virtual destructor.
    virtual ~Option() {}

    /// Evaluate and set associated option.
    virtual void Evaluate( const size_t argc, const char* argv[], size_t& index );

    /// Virtual function that returns an XML tree describing this option.
    virtual mxml_node_t* MakeXML(  mxml_node_t *const parent ) const;

    /// Return a textual description of the parameter associated with this option
    virtual std::string GetParamTypeString() const;

    /// Format additional help information (e.g., default values).
    virtual std::ostringstream& PrintHelp( std::ostringstream& fmt /*!< Stream that the additional help information is formatted into*/ ) const;

    /// Format additional help information (e.g., default values).
    virtual void PrintWiki() const;

  protected:
    /// Pointer to associated variable.
    T* Var;

    /// Pointer to (optional) flag that will be set when this option occurs.
    bool* Flag;
  };

  /// Command line option with list argument: repeated calls will add to list
  template<class T>
  class List : 
    /// Inherit from generic command line item.
    public Item
  {
  public:
    /// Constructor.
    List( std::list<T>& list ) : m_pList( &list ) {}

    /// Virtual destructor.
    virtual ~List() {}
    
    /// Evaluate and set associated option.
    virtual void Evaluate( const size_t argc, const char* argv[], size_t& index );

    /// Virtual function that returns an XML tree describing this option.
    virtual mxml_node_t* MakeXML(  mxml_node_t *const parent ) const;

  private:
    /// Pointer to associated variable.
    std::list<T>* m_pList;
  };

  /** Command line option with vector argument.
   *\note For backward compatibility, repeated use of a vector option appends
   * the subsequent vector elements onto the previously set ones. That is, the
   * vector is only cleared upon the first option use, to remove potentially
   * present default elements.
   */
  template<class T>
  class Vector : 
    /// Inherit from generic command line item.
    public Item
  {
  public:
    /// Constructor.
    Vector( std::vector<T>& vector ) : m_pVector( &vector ), m_HasBeenUsed( false ) {}
    
    /// Virtual destructor.
    virtual ~Vector() {}

    /// Evaluate and set associated option.
    virtual void Evaluate( const size_t argc, const char* argv[], size_t& index );

    /// Virtual function that returns an XML tree describing this option.
    virtual mxml_node_t* MakeXML(  mxml_node_t *const parent ) const;

    /// Return a textual description of the parameter associated with this option
    virtual std::string GetParamTypeString() const;

  private:
    /// Pointer to associated variable.
    std::vector<T>* m_pVector;

    /// Has this vector option been used already? This is so we clear the vector exactly once.
    bool m_HasBeenUsed;
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
    
    /// Virtual destructor.
    virtual ~Callback() {}

    /// Evaluate and set associated option.
    virtual void Evaluate( const size_t argc, const char* argv[], size_t& index );

    /// Virtual function that returns an XML tree describing this option.
    virtual mxml_node_t* MakeXML(  mxml_node_t *const parent ) const;

    /// Virtual function returns a string that describes the parameter type associated with this callback.
    virtual std::string GetParamTypeString() const;

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

  /// Non-option parameter.
  class NonOptionParameter :
    /// This is like a standalone string option, so inherit from that.
    public Option<const char*>
  {
  public:
    /// This class.
    typedef NonOptionParameter Self;

    /// Smart pointer.
    typedef SmartPointer<Self> SmartPtr;
    
    /// Superclass.
    typedef Option<const char*> Superclass;

    /// Constructor.
    NonOptionParameter( const char* *const var, const char* name, const std::string& comment, bool *const flag ) : Superclass( var, flag ), m_Name( name ), m_Comment( comment ) {};

    /// Evaluate and set associated variable.
    virtual void Evaluate( const size_t argc, const char* argv[], size_t& index );

    /// Returns an XML tree describing this parameter.
    virtual mxml_node_t* MakeXMLWithIndex( mxml_node_t *const parent /*!< Parent in the XML tree for the new node.*/,
					   const int index /*!< Running index [0,1,...] of this argument in the argument list.*/ ) const;
    
    /// Return a textual description of the parameter associated with this option
    virtual std::string GetParamTypeString() const;

    /// Format additional help information (e.g., default values).
    virtual std::ostringstream& PrintHelp( std::ostringstream& fmt /*!< Stream that the additional help information is formatted into*/ ) const
    {
      // by default, simply return stream unchanged.
      if ( this->Var )
	fmt << "\n[Default: " << *(this->Var) << "]";
      else
	fmt << "\n[There is no default for this parameter]";
      return fmt;
    }

    /// Format additional help information (e.g., default values).
    virtual void PrintWiki() const
    {
      // by default, simply return stream unchanged.
      if ( this->Var )
	StdOut << " '''[Default: " << *(this->Var) << "]'''\n";
      else
	StdOut << " '''[There is no default for this parameter]'''\n";
    }

    /// Name of this parameter.
    const char* m_Name;

    /// Comment (description) of this parameter.
    const std::string m_Comment;
  };

  /// Non-option parameter.
  class NonOptionParameterVector :
    /// This is like a standalone string option, so inherit from that.
    public Option< std::vector<std::string> >
  {
  public:
    /// This class.
    typedef NonOptionParameterVector Self;

    /// Smart pointer.
    typedef SmartPointer<Self> SmartPtr;
    
    /// Superclass.
    typedef Option< std::vector<std::string> > Superclass;

    /// Constructor.
    NonOptionParameterVector( std::vector<std::string> *pvec, const char* name, const std::string& comment, bool *const flag ) 
      : Superclass( pvec, flag ), 
	m_Name( name ), 
	m_Comment( comment ) {};

    /// Evaluate and set associated variable.
    virtual void Evaluate( const size_t argc, const char* argv[], size_t& index );

    /// Returns an XML tree describing this parameter.
    virtual mxml_node_t* MakeXMLWithIndex( mxml_node_t *const parent /*!< Parent in the XML tree for the new node.*/,
					   const int index /*!< Running index [0,1,...] of this argument in the argument list.*/ ) const;

    /// Return a textual description of the parameter associated with this option
    virtual std::string GetParamTypeString() const;

    /// Format additional help information (e.g., default values).
    virtual std::ostringstream& PrintHelp( std::ostringstream& fmt /*!< Stream that the additional help information is formatted into*/ ) const
    {
      if ( this->Var->size() )
	{
	fmt << "\n[Default: ( \"" << (*this->Var)[0] << "\"";
	for ( size_t i = 1; i < this->Var->size(); ++i )
	  fmt << ", \"" << (*this->Var)[i] << "\" ";
	fmt << ") ]";
	}
      else
	{
	fmt << "\n[Default: (empty)]";
	}
      return fmt;
    }
    
    /// Format additional help information (e.g., default values).
    virtual void PrintWiki() const
    {
      if ( this->Var->size() )
	{
	StdOut << "'''[Default: ( \"" << (*this->Var)[0] << "\"";
	for ( size_t i = 1; i < this->Var->size(); ++i )
	  StdOut << ", \"" << (*this->Var)[i] << "\" ";
	StdOut << ") ]'''\n";
	}
      else
	{
	StdOut << "'''[Default: (empty)]'''\n";
	}
    }
    
    /// Name of this parameter.
    const char* m_Name;

    /// Comment (description) of this parameter.
    const std::string m_Comment;
  };

public:
  /// Constructor using const inputs.
  CommandLine( const int properties = PROPS_NOXML );

  /// Destructor: spit out a warning if there are unused extra arguments on the command line.
  ~CommandLine();

  /// Forward declaration.
  class EnumGroupBase;

  /// Local class that connects command line options with their evaluators.
  class KeyToAction
  {
  public:
    /// Smart pointer to this class.
    typedef SmartPointer<KeyToAction> SmartPtr;

    /// Key for this key-action pair.
    const Key m_Key;

    /// Constructor.
    KeyToAction( const Key& key /*!< Key: long and/or short command line option for this action.*/,
		 const std::string& comment /*!< Command line help comment for this action.*/ ) :
      m_Key( key ),
      m_Comment( comment ),
      m_Properties( PROPS_XML )
    {}
    
    /// Virtual destructor.
    virtual ~KeyToAction() {};

    /// Test long key from command line and execute if match.
    virtual bool MatchAndExecute( const std::string& key /*!< Key (long option) from the command line.*/,
				  const size_t argc /*!< Total number of command line arguments.*/,
				  const char* argv[] /*!< Command line argument list.*/,
				  size_t& index /*!< Current index in command line list*/ ) = 0;

    /// Test short key from command line and execute if match.
    virtual bool MatchAndExecute( const char keyChar /*!< Key (long option) from the command line.*/,
				  const size_t argc /*!< Total number of command line arguments.*/,
				  const char* argv[] /*!< Command line argument list.*/,
				  size_t& index /*!< Current index in command line list*/ ) = 0;

    /// Set action properties.
    virtual void SetProperties( const long int properties );

    /// Get item properties.
    virtual long int GetProperties() const;
    
    /// Returns an XML tree describing this key and action.
    virtual mxml_node_t* MakeXML( mxml_node_t *const parent /*!< Parent in the XML tree for the new node.*/ ) const;

    /// Print help for this item.
    virtual void PrintHelp( const size_t globalIndent = 0 ) const = 0;
    
  protected:
    /// Format help for key part of this key/action..
    virtual void FormatHelp( std::ostringstream& fmt ) const;    

    /// Print help for this item.
    virtual void PrintWikiWithPrefix( const std::string& prefix = "" ) const;

    /// Get type info for action parameter (if any).
    virtual std::string GetActionTypeInfo() const 
    {
      return "";
    }
    
    /// Comment (description).
    std::string m_Comment;

    /// Group properties.
    long int m_Properties;

    /** Match two long options but be tolerant to hyphens, i.e., consider '-' and '_' the same.
     * This allows us to be tolerant with Slicer's requirement that there are no hyphens in
     * long options, while maintaining the ability to use them on the command line for
     * compatibility.
     *\return true is the two string match, or their only differences are hyphens vs. underlines.
     */
    virtual bool MatchLongOption( const std::string& key ) const;

    /// Give command line class full access.
    friend class CommandLine;

    /// Give enum group class full access.
    friend class CommandLine::EnumGroupBase;
  };

  /// Local class that connects command line options with single action evaluators.
  class KeyToActionSingle :
    /// Inherit from generic key-to-action class.
    public KeyToAction
  {
  public:
    /// Superclass.
    typedef KeyToAction Superclass;

    /// Smart pointer.
    typedef SmartPointer<KeyToActionSingle> SmartPtr;

    /// Constructor.
    KeyToActionSingle( const Key& key /*!< Key: long and/or short command line option for this action.*/,
		       Item::SmartPtr action /*!< The actual action (option, switch, callback, etc.)*/,
		       const std::string& comment /*!< Command line help comment for this action.*/ ) : 
      KeyToAction( key, comment ),
      m_Action( action )
    {}
    
    /// Virtual destructor.
    virtual ~KeyToActionSingle() {};

    /// Test long key from command line and execute if match.
    bool MatchAndExecute( const std::string& key /*!< Key (long option) from the command line.*/,
			  const size_t argc /*!< Total number of command line arguments.*/,
			  const char* argv[] /*!< Command line argument list.*/,
			  size_t& index /*!< Current index in command line list*/ );

    /// Test short key from command line and execute if match.
    bool MatchAndExecute( const char keyChar /*!< Key (long option) from the command line.*/,
			  const size_t argc /*!< Total number of command line arguments.*/,
			  const char* argv[] /*!< Command line argument list.*/,
			  size_t& index /*!< Current index in command line list*/ );

    /// Returns an XML tree describing this key and action.
    virtual mxml_node_t* MakeXML( mxml_node_t *const parent /*!< Parent in the XML tree for the new node.*/ ) const;

    /// Print help for this item.
    virtual void PrintHelp( const size_t globalIndent = 0 ) const;
    
    /// Print wiki help for this item.
    virtual void PrintWikiWithPrefix( const std::string& prefix = "" ) const;
    
    /// Get type info for action parameter (if any).
    virtual std::string GetActionTypeInfo() const
    {
      return this->m_Action->GetParamTypeString();
    }

    /// Action for simple key-action correspondence..
    Item::SmartPtr m_Action;
  };

  /// Base class for templated EnumGroup class.
  class EnumGroupBase :
    /// Inherit from STL list class.
    public std::list< SmartPointer<KeyToActionSingle> >
  {
  public:
    /// Smart pointer.
    typedef SmartPointer<EnumGroupBase> SmartPtr;

    /// Get key for the default value.
    std::string GetDefaultKey() const
    {
      for ( const_iterator it = this->begin(); it != this->end(); ++it )
	{
	if ( (*it)->m_Action->IsDefault() )
	  {
	  return std::string( (*it)->m_Key.m_KeyString );
	  }
	}
      return "";
    }
  };

  /// Local class that connects command line options with enum group evaluators.
  class KeyToActionEnum :
    /// Inherit from generic key-to-action class.
    public KeyToAction
  {
  public:
    /// Superclass.
    typedef KeyToAction Superclass;

    /// Smart pointer.
    typedef SmartPointer<KeyToActionEnum> SmartPtr;

    /** Constructor for enumeration parameter group.
     * An enumeration parameter group is a group of parameters that all modify the same variable
     * by setting it to different values. There is one command line parameter (i.e., key/action pair)
     * per value, plus a group parameter, which sets the variable based on the name of one of
     * the supported values.
     */
    KeyToActionEnum( const Key& key /*!< Key: long and/or short command line option for this action.*/,
		     EnumGroupBase::SmartPtr keyToActionEnum /*!< The definition of the enumeration keys and values.*/,
		     const std::string& comment /*!< Command line help comment for this action.*/ ) :
      KeyToAction( key, comment ),
      m_EnumGroup( keyToActionEnum )
    {}
    
    /// Virtual destructor.
    virtual ~KeyToActionEnum() {};

    /// Test long key from command line and execute if match.
    bool MatchAndExecute( const std::string& key /*!< Key (long option) from the command line.*/,
			  const size_t argc /*!< Total number of command line arguments.*/,
			  const char* argv[] /*!< Command line argument list.*/,
			  size_t& index /*!< Current index in command line list*/ );
    
    /// Test short key from command line and execute if match.
    bool MatchAndExecute( const char keyChar /*!< Key (long option) from the command line.*/,
			  const size_t argc /*!< Total number of command line arguments.*/,
			  const char* argv[] /*!< Command line argument list.*/,
			  size_t& index /*!< Current index in command line list*/ );

    /// Returns an XML tree describing this key and action.
    virtual mxml_node_t* MakeXML( mxml_node_t *const parent /*!< Parent in the XML tree for the new node.*/ ) const;

    /// Print help for this item.
    virtual void PrintHelp( const size_t globalIndent = 0 ) const;
    
    /// Print help for this item in Wiki markup.
    virtual void PrintWikiWithPrefix( const std::string& prefix = "" ) const;
    
  private:
    /// For enum parameter group, list of subkeys and action.
    EnumGroupBase::SmartPtr m_EnumGroup;

    /// Give enum group class full access.
    friend class EnumGroupBase;
  };

  /// Key-to-action list for enumeration parameter groups.
  template<class TDataType>
  class EnumGroup :
    /// Inherit from base class.
    public EnumGroupBase
  {
  public:
    /// Smart pointer to this class.
    typedef SmartPointer< EnumGroup<TDataType> > SmartPtr;
    
    /// Constructor.
    EnumGroup( TDataType *const variable /*!< The variable handled by this enum group.*/ ) :
      m_Variable( variable )
    {
    }

    /// Add switch to this group.
    Item::SmartPtr& 
    AddSwitch( const Key& key, const TDataType& value, const std::string& comment ) 
    {
      KeyToActionSingle::SmartPtr keyToAction( new KeyToActionSingle( key, Item::SmartPtr( new Switch<TDataType>( this->m_Variable, value ) ), comment ) );
      this->push_back( keyToAction );
      return keyToAction->m_Action;
    }

  private:
    /// Pointer to the variable that holds the enum value.
    TDataType *m_Variable;
  };

  /// Add an enum group.
  template<class TDataType>
  typename EnumGroup<TDataType>::SmartPtr AddEnum( const std::string& name, TDataType *const variable, const std::string& comment )
  {
    typename EnumGroup<TDataType>::SmartPtr enumGroup( new EnumGroup<TDataType>( variable ) );
    KeyToActionEnum::SmartPtr keyToAction( new KeyToActionEnum( Key( name ), enumGroup, comment ) );

    this->m_KeyActionList->push_back( keyToAction );
    this->m_KeyActionListComplete.push_back( keyToAction );

    return enumGroup;
  }

  /// Add switch.
  template<class T> 
  Item::SmartPtr 
  AddSwitch( const Key& key, T *const var, const T value, const std::string& comment ) 
  {
    return this->AddKeyAction( KeyToActionSingle::SmartPtr( new KeyToActionSingle( key, Item::SmartPtr( new Switch<T>( var, value ) ), comment ) ) )->m_Action;
  }
  
  /// Add option.
  template<class T> 
  Item::SmartPtr
  AddOption( const Key& key, T *const var, const std::string& comment, bool *const flag = NULL ) 
  {
    return this->AddKeyAction( KeyToActionSingle::SmartPtr( new KeyToActionSingle( key, Item::SmartPtr( new Option<T>( var, flag ) ), comment ) ) )->m_Action;
  }
  
  /// Add list option (put arguments from repeated uses into a FIFO list).
  template<class T> 
  Item::SmartPtr
  AddList( const Key& key, std::list<T>& list, const std::string& comment ) 
  {
    return this->AddKeyAction( KeyToActionSingle::SmartPtr( new KeyToActionSingle( key, Item::SmartPtr( new List<T>( list ) ), comment ) ) )->m_Action;
  }
  
  /// Add vector option (breaks argument into elements of a vector).
  template<class T> 
  Item::SmartPtr
  AddVector( const Key& key, std::vector<T>& vector, const std::string& comment ) 
  {
    return this->AddKeyAction( KeyToActionSingle::SmartPtr( new KeyToActionSingle( key, Item::SmartPtr( new Vector<T>( vector ) ), comment ) ) )->m_Action;
  }
  
  /// Add callback.
  Item::SmartPtr
  AddCallback( const Key& key, CallbackFunc func, const std::string& comment ) 
  {
    return this->AddKeyAction( KeyToActionSingle::SmartPtr( new KeyToActionSingle( key, Item::SmartPtr( new Callback( func ) ), comment ) ) )->m_Action;
  }

  /// Add callback with a single argument.
  template<class TArg>
  Item::SmartPtr
  AddCallback( const Key& key, void (*funcArg)( const TArg ), const std::string& comment ) 
  {
    return this->AddKeyAction( KeyToActionSingle::SmartPtr( new KeyToActionSingle( key, Item::SmartPtr( new Callback( funcArg ) ), comment ) ) )->m_Action;
  }
  
  /// Add callback with multiple arguments.
  Item::SmartPtr
  AddCallback( const Key& key, CallbackFuncMultiArg funcMultiArg, const std::string& comment ) 
  {
    return this->AddKeyAction( KeyToActionSingle::SmartPtr( new KeyToActionSingle( key, Item::SmartPtr( new Callback( funcMultiArg ) ), comment ) ) )->m_Action;
  }
  
  /// Add single non-option parameter.
  Item::SmartPtr
  AddParameter( const char* *const var, const char* name, const std::string& comment, bool *const flag = NULL ) 
  {
    NonOptionParameter::SmartPtr parameter( new NonOptionParameter( var, name, comment, flag ) );
    this->m_NonOptionParameterList.push_back( parameter );
    return parameter;
  }

  /// Add vector of non-option parameters.
  Item::SmartPtr
  AddParameterVector( std::vector<std::string>* pvec, const char* name, const std::string& comment, bool *const flag = NULL ) 
  {
    NonOptionParameterVector::SmartPtr vparameter( new NonOptionParameterVector( pvec, name, comment, flag ) );
    this->m_NonOptionParameterVectorList.push_back( vparameter );
    return  vparameter;
  }

  /// Type for key/action lists.
  typedef std::vector<KeyToAction::SmartPtr> KeyActionListType;

  /// Reference to current key/action list (current group).
  KeyActionListType* m_KeyActionList;
				    
  /// Reference to complete key/action list (all groups combined), which is used for key lookup.
  KeyActionListType m_KeyActionListComplete;
				    
  /// Add an item to applicable key/action lists.
  KeyToActionSingle::SmartPtr AddKeyAction( const KeyToActionSingle::SmartPtr& keyAction )
  {
    this->m_KeyActionList->push_back( keyAction );
    this->m_KeyActionListComplete.push_back( keyAction );
    return keyAction;
  }

  /// Add an item to applicable key/action lists.
  KeyToActionEnum::SmartPtr AddKeyAction( const KeyToActionEnum::SmartPtr& keyAction )
  {
    this->m_KeyActionList->push_back( keyAction );
    this->m_KeyActionListComplete.push_back( keyAction );
    return keyAction;
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

    /// Set group properties.
    virtual void SetProperties( const long int properties )
    {
      this->m_Properties = properties;
    }

    /// Get group properties.
    virtual long int GetProperties() const
    {
      return this->m_Properties;
    }
    
  private:
    /// Group properties.
    long int m_Properties;
  };

  /// Begin parameter group.
  KeyActionGroupType::SmartPtr& BeginGroup( const char* name, const char* description );

  /// End last parameter group.
  void EndGroup();

  /// Parse command line.
  bool Parse( const int argc, const char* argv[] ) throw( ExitException, Self::Exception );

  /// Help text indentation.
  static const int HelpTextIndent = 10;

  /// Print help text.
  void PrintHelp() const;

  /// Print help text.
  void PrintWiki() const;

  /// Get next parameter index.
  size_t GetNextIndex() const { return this->Index; }

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
   *\see 
   */
  void WriteXML() const;
  
private:
  /// Total number of arguments.
  size_t ArgC;

  /// Array of argument pointers.
  const char** ArgV;
  
  /// Global properties of the command line.
  long int m_Properties;

  /// Index of current argument.
  size_t Index;

  /// Type for group list.
  typedef std::vector<KeyActionGroupType::SmartPtr> KeyActionGroupListType;
  
  /// List of command line keys with associated actions.
  KeyActionGroupListType m_KeyActionGroupList;

  /// Type for no-option parameter list.
  typedef std::vector<NonOptionParameter::SmartPtr> NonOptionParameterListType;

  /// List of non-option parameters (i.e., "the rest of the command line").
  NonOptionParameterListType m_NonOptionParameterList;

  /// Type for no-option parameter vector list. These come after the scalar non-option parameters.
  typedef std::vector<NonOptionParameterVector::SmartPtr> NonOptionParameterVectorListType;

  /// List of non-option parameters (i.e., "the rest of the command line").
  NonOptionParameterVectorListType m_NonOptionParameterVectorList;

  /// Map type for program meta information.
  typedef std::map<ProgramProperties,std::string> ProgramPropertiesMapType;

  /// Program meta information.
  ProgramPropertiesMapType m_ProgramInfo;

  /// Add program info item to XML tree.
  mxml_node_t* AddProgramInfoXML( mxml_node_t *const parent /*!< Parent node for new entry in XML tree.*/,
				  const ProgramProperties key /*!< Key code for program property.*/,
				  const char* name /*!< Name of XML tag for this property.*/ ) const;

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

#include "cmtkCommandLineItem.txx"
#include "cmtkCommandLineOption.txx"
#include "cmtkCommandLineConvert.txx"
#include "cmtkCommandLineList.txx"
#include "cmtkCommandLineVector.txx"

#endif // #ifndef __cmtkCommandLine_h_included_
