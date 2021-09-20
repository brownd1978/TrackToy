//
// Implementation of class EdmlXMLEngine
//
//
#include "TrackToy/Detector/EdmlXMLEngine.hh"
#include "TrackToy/General/FileFinder.hh"

#include <iostream>
#include <sstream>
using namespace std;

#include <string.h>

// Set it to 1 when debugging

#define EdmlXMLEngine_DEBUG 0

namespace {

    bool
    nodeError( EdmlXMLEngine::Context* cxt,
               XMLNodePointer_t        node,
               const std::string&      msg )
    {
        cerr << "EdmlXMLEngine: "
             << msg
             << " [ at node: " << cxt->eng->GetNodeName( node ) << " ]" << endl;
        return false;
    }

    bool
    getNodeAttr( std::map<std::string, std::string >& theAttributes,
                 EdmlXMLEngine::Context*              cxt,
                 XMLNodePointer_t                     node )
    {
        theAttributes.clear( );

        XMLAttrPointer_t attr = cxt->eng->GetFirstAttr( node );
        while( attr ) {

            // Get and analyze the name of the currently selected attribute
            //
            const char* name = cxt->eng->GetAttrName ( attr );
            if( !name )
                return ::nodeError( cxt, node,
                                    "failed to get the attribute name" );

            if( !strlen( name ))
                return ::nodeError( cxt, node,
                                    "empty attribute name found" );

            if( theAttributes.count( name ))
                return ::nodeError( cxt, node,
                                    "the duplicate attribute: "+std::string( name )+" found" );

            // Get and analyze the value of the currently selected attribute
            //
            const char* value = cxt->eng->GetAttrValue( attr );
            if( !value )
                return ::nodeError( cxt, node,
                                    "failed to get a value for attribute: "+std::string( name ));

            // Register the attribute
            //
            // NOTE: The attribute names and values will be copied into
            //       the Standard strings during this step. Otherwise we may
            //       run into troubles when attempting to free the memory
            //       allocated to the attribute data structures.
            //
            theAttributes[name] = value;

            attr = cxt->eng->GetNextAttr( attr );
        }
        cxt->eng->FreeAllAttr( node );

        return true;
    }
}

EdmlXMLEngine::EdmlXMLEngine( )
{
#if EdmlXMLEngine_DEBUG
    cerr << this->context2string( ) << " [c-tor]" << endl;
#endif
}

EdmlXMLEngine::~EdmlXMLEngine( )
{
#if EdmlXMLEngine_DEBUG
    cerr << this->context2string( ) << " [d-tor]" << endl;
#endif
    this->reset( );
}

void
EdmlXMLEngine::reset( )
{
#if EdmlXMLEngine_DEBUG
    cerr << this->context2string( ) << " [reset]" << endl;
#endif
    while( !_cxt.empty( )) {
        delete _cxt.top( );
        _cxt.pop( );
    }
    _files.clear( );
}

bool
EdmlXMLEngine::open( const char* filename )
{
#if EdmlXMLEngine_DEBUG
    cerr << this->context2string( ) << " [open]" << endl;
#endif
    this->reset( );
    return this->doOpen( filename );
}

bool
EdmlXMLEngine::move2child( bool& noMoreNodes )
{
#if EdmlXMLEngine_DEBUG
    cerr << this->context2string( ) << " [move2child]" << endl;
#endif
    noMoreNodes = false;    // report errors by default

    if( _cxt.empty( )) {
        cerr << "EdmlXMLEngine::move2child: no parsing attempt made yet" << endl;
        return false;
    }
    Context* cxt = _cxt.top();

    // Using the ROOT node as a parent of it's the very first request
    // in the current context.
    //
    XMLNodePointer_t last = ( _cxt.top( )->last.size() > 0 ? cxt->last.top() : cxt->root );
    XMLNodePointer_t node = cxt->eng->GetChild( last );
    cxt->eng->SkipEmpty( node );
    if( node ) {

        // Check if this is <include file="name"> directive.
        //
        if( 0 == strcmp( "include", cxt->eng->GetNodeName( node ))) {

            std::map<std::string, std::string > attr;
            if( !::getNodeAttr( attr,
                                cxt,
                                node )) return false;

            if( !attr.count( "file" ))
                return ::nodeError( cxt, node,
                                    "missing attribute 'file'" );

            // Make <include> the last visited node before diving down into
            // the nested file. We may need this node to proceed with the
            // search after returning from the nested context.
            //
            cxt->last.push( node );
            if( !this->doOpen(  AppFileName(attr["file"].c_str()).pathname().c_str(), true )) {
                return false;
            }

            // Try searching down to the nested file
            //
            if( this->move2child( noMoreNodes )) return true;
            if( !noMoreNodes ) return false;

            // The nested file is completelly empty. So we need to remove
            // its context from the stack and proceed to the next (to 'include')
            // node instead.
            //
            // NOTE: This algorithm is designed to allow more than one 'include'
            //       statement to follow each other.
            //
            delete _cxt.top();
            _cxt.pop( );

            // Let's proceed horizontally to the next node (if any) instead
            //
            return this->move2next( noMoreNodes );

        } else {

            // This is what we've been looking for
            //
            cxt->last.push( node );
            return true;
        }
    }

    // We're here because we're ran out of children in the current or
    // any nested context(s).
    //
    noMoreNodes = true;

    return false;
}

bool
EdmlXMLEngine::move2next( bool& noMoreNodes )
{
#if EdmlXMLEngine_DEBUG
    cerr << this->context2string( ) << " [move2next]" << endl;
#endif
    noMoreNodes = false;
    if( !this->checkContext( "move2next" )) return false;

    Context* cxt = _cxt.top();

    XMLNodePointer_t node = cxt->eng->GetNext( cxt->last.top());
    cxt->eng->SkipEmpty( node );
    while( node ) {

        // Check if this is <include file="name"> directive. Then proceed
        // down to that file.
        //
        if( 0 == strcmp( "include", cxt->eng->GetNodeName( node ))) {

            std::map<std::string, std::string > attr;
            if( !::getNodeAttr( attr,
                                cxt,
                                node )) return false;

            if( !attr.count( "file" ))
                return ::nodeError( cxt, node,
                                    "missing attribute 'file'" );

            // Make <include> the last visited node before diving down into
            // the nested file. We may need this node to proceed with the
            // search after returning from the nested context.
            //
            cxt->last.top( ) = node;

            if( !this->doOpen( AppFileName(attr["file"].c_str()).pathname().c_str(), true )) {
                return false;
            }

            // Try searching down to the nested file
            //
            if( this->move2child( noMoreNodes )) return true;
            if( !noMoreNodes ) return false;

            // The nested file is completelly empty. So we need to remove
            // its context from the stack and proceed to the next (to 'include')
            // node instead.
            //
            // NOTE: This algorithm is designed to allow more than one 'include'
            //       statement to follow each other.
            //
            delete _cxt.top();
            _cxt.pop( );

            node = cxt->eng->GetNext( node );
            cxt->eng->SkipEmpty( node );

        } else {

            // This is what we've been looking for
            //
            cxt->last.top( ) = node;
            return true;
        }
    }
    cxt->last.pop( );

    // Check if we're in the "included" context. If so then we
    // should be handling it differently from the main one.
    //
    if( _cxt.size( ) > 1 ) {

        // If it's the root node then pop the whole context out of the stack
        // and try to see if there is a child in the previous context.
        //
        if( cxt->last.size( ) == 0 ) {
            delete _cxt.top( );
            _cxt.pop( );
            return this->move2next( noMoreNodes );
        }
    }
    noMoreNodes = true;
    return false;
}

const char*
EdmlXMLEngine::getNodeName( )
{
    if( !this->checkContext( "getNodeName" )) return 0;
    Context* cxt = _cxt.top( );
    XMLNodePointer_t node = cxt->last.top( );
    return cxt->eng->GetNodeName( node );
}

const char*
EdmlXMLEngine::getNodeContent( )
{
    if( !this->checkContext( "getNodeContent" )) return 0;
    Context* cxt = _cxt.top( );
    XMLNodePointer_t node = cxt->last.top( );
    return cxt->eng->GetNodeContent( node );
}

bool
EdmlXMLEngine::getNodeAttr( std::map<std::string, std::string >& theAttributes )
{
    if( !this->checkContext( "getNodeAttr" )) return false;

    Context* cxt = _cxt.top( );
    XMLNodePointer_t node = cxt->last.top( );

    return ::getNodeAttr( theAttributes,
                          cxt,
                          node );
}

bool
EdmlXMLEngine::doOpen( const char* filename,
                       const bool  included )
{
#if EdmlXMLEngine_DEBUG
    cerr << this->context2string( ) << " [doOpen]" << endl;
    cerr << "doOpen: filename = " << filename << ", included = " << included << endl;
#endif
    if( !filename ) {
        cerr << "EdmlXMLEngine::doOpen(): zero pointer instead of a valid file name" << endl;
        return false;
    }
    if( _files.count( filename )) {
        cerr << "EdmlXMLEngine::doOpen(): the file has already been parsed: " << filename << endl;
        return false;
    }
    TXMLEngine* xmlengine = new TXMLEngine( );
    XMLDocPointer_t xmldoc = xmlengine->ParseFile( filename );
    if( !xmldoc ) {
        cerr << "EdmlXMLEngine::doOpen(): failed to parse/validate: " << filename << endl;
        delete xmlengine;
        return false;
    }
    XMLNodePointer_t rootnode = xmlengine->DocGetRootElement( xmldoc );
    if( !rootnode ) {
        cerr << "EdmlXMLEngine::doOpen(): failed to get the root node in file: " << filename << endl;
        xmlengine->FreeDoc( xmldoc );
        delete xmlengine;
        return false;
    }
    if( 0 != strcmp( xmlengine->GetNodeName( rootnode ), "edml" )) {
        cerr << "EdmlXMLEngine::doOpen(): not an EDML document: " << filename << endl;
        xmlengine->FreeDoc( xmldoc );
        delete xmlengine;
        return false;
    }

    // For those files which are included we also want an extra tag to be
    // present placed as the only child of <edml>:
    //
    // <edml>
    //   <included>
    //     ..
    //   </included>
    // </edml>
    //
    if( included ) {

        // Make sure we have <included> as the first child of <edml>
        //
        XMLNodePointer_t node = xmlengine->GetChild( rootnode );
        xmlengine->SkipEmpty( node );
        if( !node ) {
            cerr << "EdmlXMLEngine::doOpen(): failed to get the <included> node in file: " << filename << endl;
            xmlengine->FreeDoc( xmldoc );
            delete xmlengine;
            return false;
        }
        if( 0 != strcmp( xmlengine->GetNodeName( node ), "included" )) {
            cerr << "EdmlXMLEngine::doOpen(): this EDML document can't be used as included: " << filename << endl;
            xmlengine->FreeDoc( xmldoc );
            delete xmlengine;
            return false;
        }

        // Make sure <included> is the only child of <edml>
        //
        XMLNodePointer_t next = xmlengine->GetNext( node );
        xmlengine->SkipEmpty( next );
        if( !node ) {
            cerr << "EdmlXMLEngine::doOpen(): extra tag: " << xmlengine->GetNodeName( next )
                 << " following <included> found in file: " << filename << endl;
            xmlengine->FreeDoc( xmldoc );
            delete xmlengine;
            return false;
        }

        // Switch the root node
        //
        rootnode = node;
    }
    _cxt.push( new Context( xmlengine,
                            xmldoc,
                            rootnode ));
    _files[filename] = std::string();

    return true;
}

bool
EdmlXMLEngine::checkContext( const char* theCaller ) const
{
    if( _cxt.empty( )) {
        cerr << "EdmlXMLEngine::" << theCaller << ": no parsing attempt made yet" << endl;
        return false;
    }
    if( _cxt.top( )->last.empty( )) {
        cerr << "EdmlXMLEngine::" << theCaller << ": no node selection made or parsing is over" << endl;
        return false;
    }
    return true;
}

std::string
EdmlXMLEngine::context2string( ) const
{
    ostringstream o;

    if( _cxt.empty( ))
        o << "none::";
    else {
        o << _cxt.size( ) - 1 << "::";
        if( _cxt.top( )->last.empty( ))
            o << "root::";
        else {
            o << _cxt.top( )->last.size( ) - 1 << "::"
              << _cxt.top( )->eng->GetNodeName( _cxt.top( )->last.top( ));
        }
    }
    return o.str( );
}
