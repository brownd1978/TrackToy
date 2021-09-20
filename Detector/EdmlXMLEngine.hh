/* Class: EdmlXMLEngine
 *
 * A facade to TXMLEngine to allow parsing nested EDML/XML documents
 * included from each other.
 *
 * NOTE: This is not a generic XML parser extension! The class only
 *       works with EDML documents.
 *
 * The nested files are included from any location within the body
 * of an EDML document using the following EDML directive:
 *
 *   <?xml version="1.0" encoding="UTF-8" ?>
 *   <edml>
 *     ..
 *     <include file="nested_edml_document.xml" />
 *     ..
 *   </edml>
 *
 * And here is a syntax for the nested EDML documents:
 *
 *   <?xml version="1.0" encoding="UTF-8" ?>
 *   <edml>
 *     <included>
 *       ..whatever is included..
 *     </included>
 *   </edml>
 *
 * Note the following rules for using the nested EDML files:
 *
 * - the contents of a nested file will be translated in a context of
 *   an EDML tag from which (body) the include directive was called.
 *
 * - nested files can't be parsed as top-level EDML documents. Any
 *   attempt to do so will result in a failure and error messages reported
 *   by the translator.
 *
 * - nested EDML documents may also include other nested files. There is no
 *   limit on a depth of the include process.
 *
 * - direct or indirect loops are not allowed. If the one will be detected
 *   then the translator will abort the translation and issue a error message.
 */
#ifndef EdmlXMLEngine_HH
#define EdmlXMLEngine_HH

#include "TXMLEngine.h"

#include <string>
#include <stack>
#include <map>

/// A front-end to TXMLEngine to parse nested EDML documents
/**
  */
class EdmlXMLEngine {

private:

    /* These methods and are prohibited because of the instance-level cache
     * which can't be trivially cloned in its current implementation.
     */
    EdmlXMLEngine( const EdmlXMLEngine& );
    EdmlXMLEngine& operator=( const EdmlXMLEngine& );

public:

    /// Normal c-tor

    EdmlXMLEngine( );

    /// The d-tor

    virtual ~EdmlXMLEngine( );

    /// Begin parsing
    /**
      * A top-level EDML document name is expected. This is the document which
      * will be allowed to include other (nested) EDML documents, and so forth.
      *
      * Calling this method will restart the parsing. Therefore if any parsing
      * had been done before then it will finish and all relevant data structures
      * will get released.
      *
      * Upon its successfull completion the method will return true.
      */
    bool open( const char* filename );

    /// Reset the current context
    /**
      * Calling this method would result in resetting the current context
      * as if no parsing has been done yet.
      */
    void reset( );

    /// Advance the parser's state to the first child of the current node
    /**
      * If successfull the method will change the current selection and return
      * 'true'. If no child of the current node exists or in case of any error
      * the method will return 'false', and the current selection won't change.
      * The ambiguity can be resolved by looking at a value of the method's
      * parameter 'noMoreNodes' which will be set to 'true' if no child exists,
      * or to 'false' to indicate a error. 
      */
    bool move2child( bool& noMoreNodes );

    /// Advance the parser's state to the next (to the current) node
    /**
      * If successfull the method will change the current selection and return
      * 'true'. If no next (to the current) node exists or in case of any error
      * the method will return 'false', and the current selection won't change.
      * The ambiguity can be resolved by looking at a value of the method's
      * parameter 'noMoreNodes' which will be set to 'true' if no next node exists,
      * or to 'false' to indicate a error.
      */
    bool move2next( bool& noMoreNodes );

    /// Get a name of the last selected node
    /**
      * The method will return 0 if no valid node context is available.
      */
    const char* getNodeName( );

    /// Get attributes of the last selected node
    /**
      * The method will return false if no valid node context is available.
      */
    bool getNodeAttr( std::map<std::string, std::string >& theAttributes );

    /// Get a content of the last selected node
    /**
      * The method will return 0 if no valid node context is available.
      */
    const char* getNodeContent( );

private:

    /// Make an actual file opening
    /**
      * The specified file is open and checked to ensure that it has a basic EDML
      * structure valid for the current context. The second parameter is meant to
      * tell the method how we would like to treat the file, if it's the main or
      * included EDML document.
      *
      * @return true if successful
      */
    bool doOpen( const char* filename,
                 const bool  included=false );

    /// Check if the current context is valid

    bool checkContext( const char* theCaller ) const;

    /// Return a string representation of the current context
    /**
      * The format is:
      *
      *   {<context_id>|none}::[{<last_id>:root}::[<node_name>]]
      */

    std::string context2string( ) const;

public:

    // The current context stack.
    //
    // NOTE: This class is made public to allow using it from the anonimous
    //       namespace in the implementation file of the current class.
    //
    struct Context {
        Context( TXMLEngine*      e,
                 XMLDocPointer_t  d,
                 XMLNodePointer_t r ) :
            eng (e),
            doc (d),
            root(r)
        {}
        ~Context( )
        {
            eng->FreeDoc( doc );
            delete eng;
        }
        TXMLEngine*      eng;
        XMLDocPointer_t  doc;
        XMLNodePointer_t root;

        std::stack<XMLNodePointer_t > last;
    };

private:

    std::stack<Context* > _cxt;

    // Files ever open by the parser. The key to the map is a nested
    // file, and a value is a file from which the nested file was
    // included. The top-level EDML file (the one called first)
    // won't have any value (actually it will have an empty one).
    //
    // This map is here for three reasons:
    //
    // 1. to know a context of the parsing
    // 2. to prevent cycles in the included files
    // 3. to prevent the same file to be included more than once
    //
    std::map<std::string, std::string > _files;
};

#endif // EdmlXMLEngine_HH
