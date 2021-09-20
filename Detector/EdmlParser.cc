// implementation of class EdmlParser
//

#include "TrackToy/Detector/EdmlParser.hh"
#include "TrackToy/Detector/EdmlXMLEngine.hh"
#include "TrackToy/Detector/EdmlMeasuringDevice.hh"
#include "TrackToy/Detector/EdmlCylDetElement.hh"
#include "TrackToy/Detector/EdmlRingDetElement.hh"
#include "TrackToy/Detector/EdmlDetector.hh"
#include "TrackToy/Detector/EdmlDetVolume.hh"
#include "TF1.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <sstream>

#include <map>
#include <set>

#include <assert.h>
#include <string.h>
#include <boost/tokenizer.hpp>

static const unsigned maxnpar(10);
const char* parnames[maxnpar] = { "par0", "par1", "par2", "par3", "par4", "par5", "par6", "par7", "par8", "par9"};

struct Measures {

  void add( const EdmlMeasuringDevice& e )
  {
    assert( !exists( e.name( )));
    devices.push_back( e );
    devicesByName.insert( e.name( ));
  }
  bool exists( const std::string& name ) const
  {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(",;:| ");
    tokenizer tokens(name, sep);
    bool retval(true);
    for (tokenizer::iterator titer = tokens.begin();titer != tokens.end(); ++titer){
      retval &= 0 != devicesByName.count( *titer );
    }
    return retval;
  }
  bool empty( ) const
  {
    return devicesByName.empty( );
  }
  std::vector<EdmlMeasuringDevice > devices; // to preserve an ordere of devices
  std::set<std::string > devicesByName;       // to allow quick lookup of known names
};

struct Volume {

  Volume( const std::string& vname, double dz ) :
    name(vname), zshift(dz)
  {}
  void add( EdmlDetElement* e )
  {
    assert( !exists( std::make_pair(e->name( ),e->id())));
    elements.push_back( e );
    elementsByName.insert( std::make_pair(e->name( ),e->id()));
  }
  bool exists( const std::pair<std::string,int>& ename ) const
  {
    return 0 != elementsByName.count( ename );
  }
  std::string name;
  std::vector<EdmlDetElement* > elements; // to preserve an ordere of elements
  std::set<std::pair<std::string,int> > elementsByName;  // to allow quick lookup of known names
  double zshift; // shift of this volume in Z
};

struct Detector {

  ~Detector( )
  {
    for( std::map<std::string, Volume* >::iterator itr = volumesByName.begin( );
        itr != volumesByName.end( ); ++itr ) delete itr->second;
  }
  void add( const std::string& vname, double zshift )
  {
    assert( !exists( vname ));
    volumes.push_back( vname );
    volumesByName[vname] = new Volume( vname, zshift );

  }
  void add( const std::string& vname, EdmlDetElement* e )
  {
    assert( exists( vname ));
    volumesByName[vname]->add( e );
  }
  bool exists( const std::string& vname ) const
  {
    return 0 != volumesByName.count( vname );
  }
  bool exists( const std::string& vname,
      const std::string& ename,
      const int eid ) const
  {
    assert( exists( vname ));
    std::map<std::string, Volume* >::const_iterator itr = volumesByName.find( vname );
    return itr->second->exists( std::make_pair(ename,eid) );
  }
  std::string name;
  std::vector<std::string > volumes;              // to preserve an ordere of volume names
  std::map<std::string, Volume* > volumesByName;  // volume descriptions
};

struct Config {
  void add( const std::string& path,
      const std::string& type,
      const std::string& value )
  {
    parameters[path] = std::pair<std::string,
    std::string >( type, value );
  }
  bool exists( const std::string& path ) const
  {
    return 0 != parameters.count( path );
  }
  std::map<std::string,                           // path
    std::pair<std::string,                 // type
    std::string > > parameters;  // value
};

struct Loop {
  int loopcount, loopstart, loopstep;
  std::string loopvar;
  // all the explicit lines of this loop
  std::vector<std::pair<std::string,std::map<std::string,std::string> > > lines;
  // nested loops
  std::vector<Loop> nestedloops;
};


// translate attribute values using formulas, variables, definitions, or literal values
std::string
EdmlParser::translate_attr_value( const std::string& value) const {
  static const unsigned slen(100);
  char cstr[slen];
  std::string retval;
  // start by looking for a matching formula
  std::map<std::string,Formula>::const_iterator iform = _formulas.find(value);
  if(iform != _formulas.end()){
    Formula form=iform->second;
    // set the parameters of the function according to the current variable values
    for(int ipar=0;ipar<form.parameters.size();ipar++){
      std::map<std::string,Variable>::const_iterator ivar = _variables.find(form.parameters[ipar]);
      if(ivar != _variables.end()){
        Variable var = ivar->second;
        if(var._type==Variable::intvar)
          form.formula->SetParameter(ipar,var._ival);
        else if(var._type==Variable::floatvar)
          form.formula->SetParameter(ipar,var._fval);
        else
          cerr << "illegal variable type " << var._type << " for function" << endl;
      } else
        cerr << "Can't find variable " << form.parameters[ipar] << endl;
    }
    // evalue the function and pack it in.  Note the argument is irrelevant, we vary the _parameters_
    double funval = form.formula->Eval(0.0);
    snprintf(cstr,slen,"%f",funval);
    retval = std::string(cstr);
    return retval;
  }
  // variables
  std::map<std::string,Variable>::const_iterator ivar = _variables.find(value);
  if(ivar != _variables.end()){
    Variable var =ivar->second;
    if(var._type==Variable::intvar){
      snprintf(cstr,slen,"%d",var._ival);
      retval = std::string(cstr);
    } else if(var._type==Variable::floatvar){
      snprintf(cstr,slen,"%f",var._fval);
      retval = std::string(cstr);
    } else if(var._type==Variable::stringvar){
      retval = var._sval;
    }
    return retval;
  }
  // constants
  std::map<std::string,std::string>::const_iterator iconst = _constants.find(value);
  if(iconst != _constants.end())
    return iconst->second;
  else
    // no match: return the literal
    return value;
}

bool
EdmlParser::translate_loop_attr(const std::map<std::string,std::string>& attr,Loop& loop) {
  std::map<std::string, std::string >::const_iterator itr;
  // find needed loop specification; first loop variable
  itr = attr.find( "variable" );
  if( attr.end() == itr ) {
    cerr << "error: loop variable declaration not found " << endl;
    return false;
  }
  loop.loopvar = itr->second;
  // loop count
  itr = attr.find( "count" );
  if( attr.end() == itr ) {
    cerr << "error: loop count declaration not found " << endl;
    return false;
  }
  if(1 != sscanf( translate_attr_value(itr->second).c_str(), "%d", &loop.loopcount )) {
    cerr << "error: the value: '" << itr->second << " is not an integer " << endl;
    return false;
  }
  // loop start
  itr = attr.find( "startvalue" );
  if( attr.end() == itr ) {
    cerr << "error: loop starting value declaration not found " << endl;
    return false;
  }
  if(1 != sscanf( translate_attr_value(itr->second).c_str(), "%d", &loop.loopstart )) {
    cerr << "error: the value: '" << itr->second << " is not an integer " << endl;
    return false;
  }
  // loop step (optional)
  itr = attr.find( "stepvalue" );
  if( attr.end() == itr ) {
    loop.loopstep=1;
  } else {
    if(1 != sscanf( translate_attr_value(itr->second).c_str(), "%d", &loop.loopstep )) {
      cerr << "error: the value: '" << itr->second << " is not an integer " << endl;
      return false;
    }
  }
  return true;
}


// Helper functions for translating attributes of various types
//
bool EdmlParser::translate_attr_float( float&                                     value,
    const std::map<std::string, std::string >& attr,
    const char*                                name )
{
  std::map<std::string, std::string >::const_iterator itr = attr.find( name );
  if( attr.end() == itr ) {
    //            cerr << "error: missing attribute: '" << name << "'" << endl;
    return false;
  }
  if( 1 != sscanf( translate_attr_value(itr->second).c_str(), "%f", &value )) {
    cerr << "error: the value: '" << itr->second
      << "' of attribute: '" << name << "'"
      << " is not of 'float' type." << endl;
    return false;
  }
  return true;
}

bool EdmlParser::translate_attr_float( EdmlAnyTypeDict<std::string>&              dict,
    const std::map<std::string, std::string >& attr,
    const char*                                name )
{
  float value = 0.0;
  if( translate_attr_float( value, attr, name ))
    return dict.insert( name, value );
  return false;
}

bool EdmlParser::translate_attr_float_if_avail( float&                                     value,
    const std::map<std::string, std::string >& attr,
    const char*                                name )
{
  std::map<std::string, std::string >::const_iterator itr = attr.find( name );
  if( attr.end() != itr ) {
    if( 1 != sscanf( translate_attr_value(itr->second).c_str(), "%f", &value )) {
      cerr << "error: the value: '" << itr->second
        << "' of attribute: '" << name << "'"
        << " is not of 'float' type." << endl;
      return false;
    }
  }
  return true;
}

bool EdmlParser::translate_attr_int( int&                                     value,
    const std::map<std::string, std::string >& attr,
    const char*                                name )
{
  std::map<std::string, std::string >::const_iterator itr = attr.find( name );
  if( attr.end() == itr ) {
    //            cerr << "error: missing attribute: '" << name << "'" << endl;
    return false;
  }
  if( 1 != sscanf( translate_attr_value(itr->second).c_str(), "%d", &value )) {
    cerr << "error: the value: '" << itr->second
      << "' of attribute: '" << name << "'"
      << " is not of 'int' type." << endl;
    return false;
  }
  return true;
}

bool EdmlParser::translate_attr_int_if_avail( int&                                     value,
    const std::map<std::string, std::string >& attr,
    const char*                                name )
{
  std::map<std::string, std::string >::const_iterator itr = attr.find( name );
  if( attr.end() != itr ) {
    if( 1 != sscanf( translate_attr_value(itr->second).c_str(), "%d", &value )) {
      cerr << "error: the value: '" << itr->second
        << "' of attribute: '" << name << "'"
        << " is not of 'int' type." << endl;
      return false;
    }
  }
  return true;
}

bool EdmlParser::translate_attr_float_if_avail( EdmlAnyTypeDict<std::string>&              dict,
    const std::map<std::string, std::string >& attr,
    const char*                                name )
{
  float value = 0.0;
  if( translate_attr_float_if_avail( value, attr, name ))
    return dict.insert( name, value );
  return false;
}

bool EdmlParser::translate_attr_string( std::string&                               value,
    const std::map<std::string, std::string >& attr,
    const char*                                name )
{
  std::map<std::string, std::string >::const_iterator itr = attr.find( name );
  if( attr.end() == itr ) {
    //            cerr << "error: missing attribute: '" << name << "'" << endl;
    return false;
  }
  value = translate_attr_value(itr->second).c_str();
  return true;
}

bool EdmlParser::translate_attr_string( EdmlAnyTypeDict<std::string>&              dict,
    const std::map<std::string, std::string >& attr,
    const char*                                name )
{
  std::string value = "";
  if( translate_attr_string( value, attr, name ))
    return dict.insert( name, value );
  return false;
}

bool EdmlParser::translate_attr_string_if_avail( std::string&                               value,
    const std::map<std::string, std::string >& attr,
    const char*                                name )
{
  std::map<std::string, std::string >::const_iterator itr = attr.find( name );
  if( attr.end() != itr ){
    value = translate_attr_value(itr->second);
  }
  return true;
}

bool EdmlParser::translate_attr_string_if_avail( EdmlAnyTypeDict<std::string>&              dict,
    const std::map<std::string, std::string >& attr,
    const char*                                name )
{
  std::string value = "";
  if( translate_attr_string_if_avail( value, attr, name ))
    return dict.insert( name, value );
  return false;
}

// Helper containers to store intermediate results of the translation.
//
// OWNERSHIP NOTE: The containers do NOT own the object stored
//                 by pointers! These objects are going to be
//                 extracted and used by the parser.
//


// Translate a section with descriptions of measuring devices. Put results
// into the container (the first parameter of the function).
//
// NOTE: This method can be called more than once. In that case it will
//       simply accumulate results in the output container.
//
bool EdmlParser::translate_measures( Measures&    measures,
    EdmlXMLEngine& xmlengine )
{
  bool noMoreNodes = false;

  if( xmlengine.move2child( noMoreNodes )) {

    do {

      const char* node = xmlengine.getNodeName( );
      if( !node ) return false;

      if( 0 != strcmp( node, "device" )) {
        cerr << "EdmlParser: unknown tag <" << node << "> in section <measures>" << endl;
        return false;
      }

      // Get attributes of the current device
      //
      std::map<std::string, std::string > attr;
      if( !xmlengine.getNodeAttr( attr )) return false;

      // Attributes of the base type
      //
      std::string type;
      if( !translate_attr_string( type, attr, "type" )) {
        cerr << "EdmlParser: missing attribute for a type of a measuring device" << endl;
        return false;
      }
      std::string name;
      if( !translate_attr_string( name, attr, "name" )) {
        cerr << "EdmlParser: missing name of a measuring device" << endl;
        return false;
      }
      // default time window of 1e-6
      float stwin(1e-6);
      translate_attr_float_if_avail( stwin, attr, "sensitiveTimeWindow");


      if( measures.exists( name )) {
        cerr << "EdmlParser: measuring device name: " << name
          << " was already defined" << endl;
        return false;
      }

      // Build and register an object of the found type
      //
      EdmlMeasuringDevice device( type, name, stwin );

      if( "DoubleSideSiStrips" == type ) {

        if(
            !translate_attr_float( device, attr, "phiResModel"   ) ||
            !translate_attr_float( device, attr, "phiResPar0"   ) ||
            !translate_attr_float( device, attr, "phiPTRatio"   ) ||
            !translate_attr_float( device, attr, "phingFrac"   ) ||
            !translate_attr_float( device, attr, "phingFactor"   ) ||
            !translate_attr_float( device, attr, "phiEff" ) ||
            !translate_attr_float( device, attr, "zResModel"   ) ||
            !translate_attr_float( device, attr, "zResPar0"   ) ||
            !translate_attr_float( device, attr, "zPTRatio"   ) ||
            !translate_attr_float( device, attr, "zngFrac"   ) ||
            !translate_attr_float( device, attr, "zngFactor"   ) ||
            !translate_attr_float( device, attr, "zEff" ) ||
            !translate_attr_float( device, attr, "NormRes"  )      ) {
          cerr << "EdmlParser: failed to translate attributes of measuring device "
            << " of type: " << type << " and name: " << name << endl;
          return false;
        }

      } else if( "SiPixels" == type ) {

        if(
            !translate_attr_float( device, attr, "phiPTRatio"   ) ||
            !translate_attr_float( device, attr, "zPTRatio"   ) ||
            !translate_attr_float( device, attr, "phiResModel"   ) ||
            !translate_attr_float( device, attr, "phiResPar0"   ) ||
            !translate_attr_float_if_avail( device, attr, "phiResPar1"   ) ||
            !translate_attr_float_if_avail( device, attr, "phiResPar2"   ) ||
            !translate_attr_float_if_avail( device, attr, "phiResPar3"   ) ||
            !translate_attr_float_if_avail( device, attr, "phiResPar4"   ) ||
            !translate_attr_float_if_avail( device, attr, "phiResPar5"   ) ||
            !translate_attr_float_if_avail( device, attr, "phiResPar6"   ) ||
            !translate_attr_float( device, attr, "zResModel"   ) ||
            !translate_attr_float( device, attr, "zResPar0"   ) ||
            !translate_attr_float_if_avail( device, attr, "zResPar1"   ) ||
            !translate_attr_float_if_avail( device, attr, "zResPar2"   ) ||
            !translate_attr_float_if_avail( device, attr, "zResPar3"   ) ||
            !translate_attr_float_if_avail( device, attr, "zResPar4"   ) ||
            !translate_attr_float_if_avail( device, attr, "zResPar5"   ) ||
            !translate_attr_float_if_avail( device, attr, "zResPar6"   ) ||
            !translate_attr_float( device, attr, "NormRes"   ) ||
            !translate_attr_float( device, attr, "ngFrac"   ) ||
            !translate_attr_float( device, attr, "ngFactor"   ) ||
            !translate_attr_float( device, attr, "Eff" )) {
              cerr << "EdmlParser: failed to translate attributes of measuring device "
                << " of type: " << type << " and name: " << name << endl;
              return false;
            }

      } else if( "DriftChamber" == type ) {

        if(
            !translate_attr_float( device, attr, "eff_par0"   ) ||
            !translate_attr_float( device, attr, "eff_par1"   ) ||
            !translate_attr_float( device, attr, "rms_par0"   ) ||
            !translate_attr_float( device, attr, "rms_par1"   ) ||
            !translate_attr_float( device, attr, "rms_par2"   ) ||
            !translate_attr_float( device, attr, "rms_par3"   ) ||
            !translate_attr_float( device, attr, "rms_par4"   ) ||
            !translate_attr_float( device, attr, "rms_par5"   ) ||
            !translate_attr_float_if_avail( device, attr, "ngFrac"   ) ||
            !translate_attr_float_if_avail( device, attr, "ngFactor"   ) ||
            !translate_attr_float( device, attr, "cell_size"   ) ||
            !translate_attr_float( device, attr, "angle" )) {
          cerr << "EdmlParser: failed to translate attributes of measuring device "
            << " of type: " << type << " and name: " << name << endl;
          return false;
        }
      } else if( "dEdx" == type ) {

        if(!translate_attr_float( device, attr, "dedx_par1" ) ||
            !translate_attr_float( device, attr, "HitType" ) ||
            !translate_attr_float( device, attr, "dedx_par2" ) ||
            !translate_attr_float( device, attr, "dedx_par3" ) ||
            !translate_attr_float( device, attr, "trunc_frac" ) ||
            !translate_attr_float_if_avail( device, attr, "svt_thresh" ) ||
            !translate_attr_float_if_avail( device, attr, "svt_fit_p0" ) ||
            !translate_attr_float_if_avail( device, attr, "svt_fit_p1" ) ||
            !translate_attr_float_if_avail( device, attr, "svt_fit_p2" ) ||
            !translate_attr_float_if_avail( device, attr, "svt_fit_p3" ) ||
            !translate_attr_float_if_avail( device, attr, "svt_fit_p4" ) ) {

          cerr << "EdmlParser: failed to translate attributes of measuring device "
            << " of type: " << type << " and name: " << name << endl;
          return false;
        }

      } else if( "DircCylinder" == type ) {
      } else if( "EmcCluster" == type ) {
      } else if( "IfrCylinder" == type ) {
      } else if( "ForwardPidDetector" == type ) {
      } else if( "FarichDetector" == type ) {
      } else {
        cerr << "EdmlParser: unknown measuring device type: " << type << endl;
        return false;
      }
      measures.add( device );


    } while( xmlengine.move2next( noMoreNodes ));
  }
  return noMoreNodes;
}

// Translate a volume with descriptions of sub-detector elements. Put results
// into the container (the first parameter of the function). The traslation
// process will also rely on the contents of a container with description
// of measuring devices which is passes as the second parameter to
// the function.
//
// NOTE: This method can be called more than once. In that case it will
//       simply accumulate results in the output container. Just keep in
//       mind that all measuring devices required for the sub-detector
//       elements discovered by the method should already be known.
//
bool EdmlParser::translate_volume( Detector&        detector,
    Measures&        measures,
    EdmlXMLEngine&     xmlengine,
    const std::string& vname )
{
  bool noMoreNodes = false;

  if( xmlengine.move2child( noMoreNodes )) {
    do {
      // The type of the element is determined by the tag's name
      const char* type = xmlengine.getNodeName( );
      if( !type ) return false;
      std::map<std::string, std::string > attr;
      if( !xmlengine.getNodeAttr( attr )) return false;
      // check for loops
      if( 0 == strcmp( type, "loop" )) {
        // first translate: this parses all the nested loops as well.
        Loop loop;
        bool tloop = translate_loop(xmlengine,attr,loop);
        if(!tloop){
          cerr << "failed translating loop" << endl;
          return tloop;
        }
        // now execute the loop
        bool xloop = execute_loop(loop,detector,measures,vname);
        if(!xloop){
          cerr << "failed executing loop" << endl;
          return xloop;
        }
      } else {
        // process the node as a detector element
        bool goodelem = translate_elem(detector,measures,vname,type,attr);
        if(!goodelem){
          cerr << "failed translating elem" << endl;
          return goodelem;
        }
      }
    } while( xmlengine.move2next( noMoreNodes ));
  }
  return noMoreNodes;
}


bool
EdmlParser::translate_elem( Detector&        detector,
    Measures&        measures,
    const std::string& vname, const char* type,const std::map<std::string,std::string>& attr){
  std::string name;
  if( !translate_attr_string( name, attr, "name" )) {
    cerr << "EdmlParser: missing name of the detector element" << endl;
    return false;
  }
  int id = -1;
  translate_attr_int_if_avail ( id, attr, "id" );

  if( detector.exists( vname, name, id )) {
    cerr << "EdmlParser: detector element: " << name << " : " << id
      << " was already defined for volume: " + vname << endl;
    return false;
  }

  // Build and register an object of the found type
  //
  EdmlDetElement* element = 0;
  if( 0 == strcmp( type, "cyl" )) {

    float zmin   = 0.0;
    float zmax   = 0.0;
    float radius = 0.0;
    float thick  = 0.0;

    std::string mat  = "";
    std::string meas = "";
    std::string type="Cylinder";

    float gap     = 0.0;
    float overlap = 0.0;

    if( !translate_attr_float          ( zmin,    attr, "zmin"    ) ||
        !translate_attr_float          ( zmax,    attr, "zmax"    ) ||
        !translate_attr_float          ( radius,  attr, "radius"  ) ||
        !translate_attr_float          ( thick,   attr, "thick"   ) ||
        !translate_attr_string_if_avail( mat,     attr, "mat"     ) ||
        !translate_attr_string_if_avail( meas,    attr, "meas"    ) ||
        !translate_attr_float_if_avail ( gap,     attr, "gap"     ) ||
        !translate_attr_float_if_avail ( overlap, attr, "overlap" ) ||
        !translate_attr_string_if_avail ( type, attr, "type" )) {


      cerr << "error: failed to translate attributes of detector element "
        << " of type: " << type << " and name: " << name << endl;
      return false;
    }

    // Make sure the measuring device is known
    //
    if( !meas.empty( )) {
      if( !measures.exists( meas )) {
        cerr << "EdmlParser: unknown measuring device: " << meas
          << " requested for detector element of type: " << type
          << " and name: " << name << endl;
        return false;
      }
    }
    element = new EdmlCylDetElement( type,
        name,
        mat,
        meas,
        id,
        thick,
        gap,
        overlap,
        zmin,
        zmax,
        radius);

  } else if( 0 == strcmp( type, "ring" )) {

    float z         = 0.0;
    float lowradius = 0.0;
    float hiradius  = 0.0;
    float thick     = 0.0;

    std::string mat  = "";
    std::string meas = "";
    std::string type = "Ring";

    float gap     = 0.0;
    float overlap = 0.0;

    if( !translate_attr_float          ( z,         attr, "z"         ) ||
        !translate_attr_float          ( lowradius, attr, "lowradius" ) ||
        !translate_attr_float          ( hiradius,  attr, "hiradius"  ) ||
        !translate_attr_float          ( thick,     attr, "thick"     ) ||
        !translate_attr_string_if_avail( mat,       attr, "mat"       ) ||
        !translate_attr_string_if_avail( meas,      attr, "meas"      ) ||
        !translate_attr_float_if_avail ( gap,       attr, "gap"     ) ||
        !translate_attr_float_if_avail ( overlap,   attr, "overlap" ) ||
        !translate_attr_string_if_avail ( type,   attr, "type" )) {


      cerr << "EdmlParser: failed to translate attributes of detector element "
        << " of type: " << type << " and name: " << name << endl;
      return false;
    }

    // Make sure the measuring device is known
    //
    if( !meas.empty( )) {
      if( !measures.exists( meas )) {
        cerr << "EdmlParser: unknown measuring device: " << meas
          << " requested for detector element of type: " << type
          << " and name: " << name << endl;
        return false;
      }
    }
    element = new EdmlRingDetElement( type,
        name,
        mat,
        meas,
        id,
        thick,
        gap,
        overlap,
        z,
        lowradius,
        hiradius);
  } else if( 0 == strcmp( type, "rect" )) {
    int orient      = -1;
    float rho       = 0.0;
    float z       = 0.0;
    float phi       = 0.0;
    float eta       = 0.0;
    float usize     = 0.0;
    float vsize     = 0.0;
    float thick     = 0.0;

    std::string mat  = "";
    std::string meas = "";
    std::string type = "Rect";

    float gap     = 0.0;
    float overlap = 0.0;

    if( !translate_attr_int ( orient, attr, "orientation" ) ||
        !translate_attr_float          ( rho,         attr, "centerrho"         ) ||
        !translate_attr_float          ( z,         attr, "centerz"         ) ||
        !translate_attr_float          ( phi, attr, "centerphi" ) ||
        !translate_attr_float          ( eta,  attr, "axisphi"  ) ||
        !translate_attr_float          ( usize, attr, "usize" ) ||
        !translate_attr_float          ( vsize,  attr, "vsize"  ) ||
        !translate_attr_float          ( thick,     attr, "thick"     ) ||
        !translate_attr_string_if_avail( mat,       attr, "mat"       ) ||
        !translate_attr_string_if_avail( meas,      attr, "meas"      ) ||
        !translate_attr_float_if_avail ( gap,       attr, "gap"     ) ||
        !translate_attr_float_if_avail ( overlap,   attr, "overlap" ) ||
        !translate_attr_string_if_avail ( type,   attr, "type" )) {


      cerr << "EdmlParser: failed to translate attributes of detector element "
        << " of type: " << type << " and name: " << name << endl;
      return false;
    }

    // Make sure the measuring device is known
    //
    if( !meas.empty( )) {
      if( !measures.exists( meas )) {
        cerr << "EdmlParser: unknown measuring device: " << meas
          << " requested for detector element of type: " << type
          << " and name: " << name << endl;
        return false;
      }
    }
    element = new EdmlRectDetElement( type,
        name,
        mat,
        meas,
        id,
        thick,
        gap,
        overlap,
        (EdmlRectDetElement::porient)orient,
        rho,z,phi,eta,usize,vsize);

  } else if( 0 == strcmp( type, "cone" )) {

    float rmin   = 0.0;
    float rmax   = 0.0;
    float tantheta = 0.0;
    float zvertex = 0.0;
    float thick  = 0.0;

    float rho1,rho2,z1,z2;

    std::string mat  = "";
    std::string meas = "";
    std::string type = "Cone";

    float gap     = 0.0;
    float overlap = 0.0;
    // we allow 2 cone representations

    bool cone1 =
      translate_attr_float( rmin,    attr, "rmin"    ) &&
      translate_attr_float( rmax,    attr, "rmax"    ) &&
      translate_attr_float( tantheta,  attr, "tantheta"  ) &&
      translate_attr_float( zvertex,  attr, "zvertex"  );

    bool cone2 =
      translate_attr_float( rho1,    attr, "rho1"    ) &&
      translate_attr_float( rho2,    attr, "rho2"    ) &&
      translate_attr_float( z1,  attr, "z1"  ) &&
      translate_attr_float( z2,  attr, "z2"  );

    if( !(cone1 || cone2) ||
        !translate_attr_float          ( thick,   attr, "thick"   ) ||
        !translate_attr_string_if_avail( mat,     attr, "mat"     ) ||
        !translate_attr_string_if_avail( meas,    attr, "meas"    ) ||
        !translate_attr_float_if_avail ( gap,     attr, "gap"     ) ||
        !translate_attr_float_if_avail ( overlap, attr, "overlap" ) ||
        !translate_attr_string_if_avail ( type, attr, "type" )) {

      cerr << "error: failed to translate attributes of detector element "
        << " of type: " << type << " and name: " << name << endl;
      return false;
    }
    // Make sure the measuring device is known
    //
    if( !meas.empty( )) {
      if( !measures.exists( meas )) {
        cerr << "EdmlParser: unknown measuring device: " << meas
          << " requested for detector element of type: " << type
          << " and name: " << name << endl;
        return false;
      }
    }
    if(cone1) {
      element = new EdmlConeDetElement( type,
          name,
          mat,
          meas,
          id,
          thick,
          gap,
          overlap,
          rmin,
          rmax,
          tantheta,
          zvertex);
    } else {
      std::pair<float,float> cyl1 = std::make_pair(rho1,z1);
      std::pair<float,float> cyl2 = std::make_pair(rho2,z2);

      element = new EdmlConeDetElement( type,
          name,
          mat,
          meas,
          id,
          thick,
          gap,
          overlap,
          cyl1,
          cyl2);
    }
  } else {
    cerr << "EdmlParser: unknown sub-detector type: " << type << endl;
    return false;
  }
  detector.add( vname, element );
  return true;
}

// Translate a section with dvolumes.
//
// NOTE: This method can be called more than once. In that case it will
//       simply accumulate results in the output container. Just keep in
//       mind that all measuring devices required for the sub-detector
//       elements discovered by the method should already be known.
//
bool EdmlParser::translate_detector( Detector&    detector,
    Measures&    measures,
    EdmlXMLEngine& xmlengine )
{
  bool noMoreNodes = false;

  if( xmlengine.move2child( noMoreNodes )) {

    do {

      const char* node = xmlengine.getNodeName( );
      if( !node ) return false;

      if( 0 != strcmp( node, "volume" )) {
        cerr << "EdmlParser: the <volume> section expected" << endl;
        return false;
      }

      // Get and analyze attributes of the volume
      //
      std::map<std::string, std::string > attr;
      if( !xmlengine.getNodeAttr( attr )) return false;

      std::string vname;
      if( !translate_attr_string( vname, attr, "name" )) {
        cerr << "EdmlParser: missing name of the volume" << endl;
        return false;
      }
      // volume parameters: for now, just the Z shift
      float zshift(0.0);
      translate_attr_float_if_avail( zshift, attr, "zshift");

      if( detector.exists( vname )) {
        cerr << "EdmlParser: volume: " << vname << " was already defined" << endl;
        return false;
      }

      // Register the new volume and translate its contents
      //
      detector.add( vname, zshift );

      if( !translate_volume( detector,
            measures,
            xmlengine,
            vname )) return false;

    } while( xmlengine.move2next( noMoreNodes ));
  }
  return noMoreNodes;
}

// Translate user-defined configuration parameters.
//
// NOTE: This method can be called more than once. In that case it will
//       simply accumulate results in the output container. The function
//       can also be called recursivelly for nested sections.
//
bool EdmlParser::translate_config( Config&          config,
    EdmlXMLEngine&     xmlengine,
    const std::string& basepath )
{
  // Stop recursion if no child sections or parameters found
  //
  bool noMoreNodes = false;

  if( xmlengine.move2child( noMoreNodes )) {

    do {

      const char* node = xmlengine.getNodeName( );
      if( !node ) return false;

      // Get and analyze common attributes of the node
      //
      std::map<std::string, std::string > attr;
      if( !xmlengine.getNodeAttr( attr )) return false;

      std::string name;
      if( !translate_attr_string( name, attr, "name" )) {
        cerr << "EdmlParser: missing name of the config element at path: " << basepath << endl;
        return false;
      }
      const std::string path = basepath + name;
      if( config.exists( path )) {
        //          cerr << "EdmlParser: config secton or parameter: " << path << " already defined" << endl;
        //          return false;
      }

      // Analyze a tag name of the current node
      //
      if( 0 == strcmp( node, "param" )) {

        // Get a type and a value (both should exist) of the parameter
        //
        std::string type;
        if( !translate_attr_string( type, attr, "type" )) {
          cerr << "EdmlParser: missing type attribute of config parameter: " << path << endl;
          return false;
        }
        config.add( path,
            type,
            xmlengine.getNodeContent( ));

      } else if( 0 == strcmp( node, "sect" )) {

        // Proceed down recursivelly to the sub-section
        //
        if( !translate_config( config,
              xmlengine,
              path + "." )) return false;

      } else {
        cerr << "EdmlParser: unknown tag <" << node << "> in config path: " << basepath << endl;
        return false;
      }

    } while( xmlengine.move2next( noMoreNodes ));
  }
  return noMoreNodes;
}

bool EdmlParser::translate( Detector& detector,
    Measures& measures,
    Config&   config,
    const char* filename )
{
  EdmlXMLEngine xmlengine;

  if( !xmlengine.open( filename )) {
      cerr << "EdmlParser: failed to parse/validate: " << filename << endl;
      return false;
    }

        // This flag is used during the subseqeunt navigation in the EDML tree.
        // Its value is set by calls to the corresponding methods of the parsing
        // engine.
        //
    bool noMoreNodes = false;    // initial value makes no sense

        // Begin the parsing by getting the very first node
        //
    if( !xmlengine.move2child( noMoreNodes )) return noMoreNodes;

    do {

      const char* node = xmlengine.getNodeName( );
      if( !node ) return false;
// strip out 'included'
      if( 0 == strcmp( node, "included" )) {
        if( !xmlengine.move2child( noMoreNodes )) return noMoreNodes;
        node = xmlengine.getNodeName( );
      }

      if( 0 == strcmp( node, "define" )) {

        if( !translate_definitions( xmlengine )) return false;

      } else if( 0 == strcmp( node, "measures" )) {

                // A description of known measuring devices and their properties.
                //
        if( !translate_measures( measures,
            xmlengine )) return false;

      } else if( 0 == strcmp( node, "detector" )) {

                // A section with a description of the detector assembly.
                //
                // NOTE: This section should be defined _After_ the one
                //       with the measuring devices. We're not cheching
                //       in the current implementation though.
                //
        std::map<std::string, std::string > attr;
        if( !xmlengine.getNodeAttr( attr )) return false;

        std::string detname;
        translate_attr_string_if_avail( detname,
          attr,
          "name" );
        if(detname.empty() && detector.name.empty()){
          cerr << "EdmlParser: missing attribute for detector name" << endl;
          return false;
        } else if (!detector.name.empty() && !detname.empty()  && detector.name != detname){
          cerr << "EdmlParser: inconsistent detector names " << detname << " " << detector.name << endl;
          return false;
        } else if(detector.name.empty()){
          detector.name = detname;
        }

        if( !translate_detector( detector,
          measures,
          xmlengine )) return false;

      } else  if( 0 == strcmp( node, "config" )) {

                // Check for the optional <config> section with user-defined
                // sections and parameters of the application.
                //
        if( !translate_config( config,
          xmlengine )) return false;
      } else {
        cerr << "EdmlParser: unknown section: " << node << endl;
        return false;
      }

    } while( xmlengine.move2next( noMoreNodes ));

    return noMoreNodes;
  }

EdmlParser::EdmlParser( const bool verbose ) :
_verbose(verbose),
  _detector(0)
  { }

EdmlParser::~EdmlParser( )
{
  this->reset( );
}

bool
EdmlParser::parse( const char* filename )
{
  if( !filename ) {
    if( _verbose )
      cout << "EdmlParser: file name expected" << endl;
    return false;
  }

    // Parse the file and (if successfull) get results into
    // the temporary containers.
    //
  Detector detector;
  Measures measures;
  Config   config;

  if( !translate( detector,
    measures,
    config,
    filename )) return false;

    // Adopt the results into the current object's context
    //
  _filenames.push_back( filename );

  if( detector.volumes.size() > 0) {
//    delete _detector;
    std::vector<EdmlDetVolume* > volumes;
    for( std::vector<std::string >::const_iterator itr = detector.volumes.begin();
    itr != detector.volumes.end(); ++itr ) {
      Volume* volume = detector.volumesByName[*itr];
      volumes.push_back( new EdmlDetVolume( volume->name, volume->elements, volume->zshift ));
    }
    if(_detector == 0){
      _detector = new EdmlDetector( detector.name, volumes );
      _devices  = measures.devices;
    } else {
      _detector->append(volumes);
      _devices.insert(_devices.end(),measures.devices.begin(),measures.devices.end());
    }
  }
  for( std::map<std::string, std::pair<std::string, std::string > >::const_iterator itr = config.parameters.begin( );
  itr != config.parameters.end( ); ++itr ) {
    _config[itr->first] = itr->second;
  }
  return true;
}

void
  EdmlParser::filenames( std::vector<std::string >& theNames ) const
{
  theNames.clear( );
  for( std::vector<std::string >::const_iterator itr = _filenames.begin( );
  itr != _filenames.end( ); ++itr ) theNames.push_back( *itr );
}

void
  EdmlParser::devices( std::vector<EdmlMeasuringDevice >& theDevices ) const
{
  theDevices = _devices;
}

void
  EdmlParser::config( std::map<std::string, std::pair<std::string, std::string > >& theParameters ) const
{
  theParameters = _config;
}

void
  EdmlParser::reset( )
{
  _filenames.clear( );
  if( _detector ) delete _detector;
  _detector = 0;
  _devices.clear( );
  _config.clear( );
}

bool
  EdmlParser::dump( std::ostream& str ) const
{
  str << "\n"
    << "[ EDML FILES TRANSLATED ]\n"
    << "\n";
  for( std::vector<std::string >::const_iterator itr = _filenames.begin( );
  itr != _filenames.end( ); ++itr ) str << "  " << *itr << "\n";

  str << "\n"
    << "[ DEFINITIONS ]\n"
    << "\n";
  for( std::map<std::string,std::string >::const_iterator itr = _constants.begin( );
  itr != _constants.end( ); ++itr ) str << " Name: " << itr->first << " Value: " << itr->second << "\n";

  str << "\n"
    << "[ MEASURING DEVICES ]\n"
    << "\n";
  for( std::vector<EdmlMeasuringDevice >::const_iterator itr = _devices.begin( );
  itr != _devices.end( ); ++itr ) str << "  " << itr->toString( ) << "\n";

  if( _detector ) {

    str << "\n"
      << "[ DETECTOR: '" << _detector->name() << "' ]\n";

    std::vector<const EdmlDetVolume* > volumes;
    _detector->volumes( volumes );

    for( std::vector<const EdmlDetVolume* >::const_iterator itr = volumes.begin( );
    itr != volumes.end( ); ++itr ) {

      const EdmlDetVolume* volume = (*itr);

      str << "\n"
        << "  [ VOLUME: '" << volume->name( ) << " Z shift " << volume->zShift() << "' ]\n"
        << "\n";
      std::vector<const EdmlDetElement* > elements;
      volume->elements( elements );

      for( std::vector<const EdmlDetElement* >::const_iterator itr = elements.begin( );
      itr != elements.end( ); ++itr ) str << "    " << (*itr)->toString( ) << "\n";
    }
  }

  str << "\n"
    << "[ CONFIGURATION PARAMETERS ]\n"
    << "\n";
  for( std::map<std::string, std::pair<std::string, std::string > >::const_iterator itr = _config.begin( );
  itr != _config.end( ); ++itr ) str << "  " << itr->first << " [" << itr->second.first << "] '" << itr->second.second << "'\n";

  str << "\n"
    << "[ END ]\n"
    << endl;

  return true;
}


bool
EdmlParser::translate_definitions(  EdmlXMLEngine& xmlengine )
{
  bool noMoreNodes = false;

  if( xmlengine.move2child( noMoreNodes )) {

    do {

      const char* node = xmlengine.getNodeName( );
      if( !node ) return false;
      std::string name;
      std::map<std::string, std::string > attr;
      if( !xmlengine.getNodeAttr( attr )) return false;

      if( 0 == strcmp( node, "constant" )) {

        if( !translate_attr_string( name, attr, "name" )) {
          cerr << "EdmlParser: missing name of the constant" << endl;
          return false;
        }

      // constant parameters: for now, just one - its value (as a string)
      //
        std::string value;
        if( !translate_attr_string( value, attr, "value" )) {
          cerr << "EdmlParser: missing value of the constant" << endl;
          return false;
        }
      // Register the new constant or replace an existing one if it was
      // defined before
      //
        _constants[name] = value;
      } else if ( 0 == strcmp( node, "variable" )){
        if( !translate_attr_string( name, attr, "name" )) {
            cerr << "EdmlParser: missing name of the variable" << endl;
            return false;
        }
        std::string type;
        if( !translate_attr_string( type, attr, "type" )) {
          cerr << "EdmlParser: missing attribute for a type of a variable" << endl;
          return false;
        }
        Variable var;
        if(type=="int")
          var._type=Variable::intvar;
        else if(type=="float")
          var._type=Variable::floatvar;
        else if(type=="string")
          var._type=Variable::stringvar;
        else{
          cerr << "Variable " << name << " has unknown type " << type << endl;
          return false;
        }
        _variables[name] = var;
        std::map<std::string,Variable>::const_iterator ifind = _variables.find(name);
        if(ifind == _variables.end())
          cerr << "error finding variable " << endl;
      } else if(0 == strcmp( node, "formula" )){
        if( !translate_attr_string( name, attr, "name" )) {
            cerr << "EdmlParser: missing name of the Formula" << endl;
            return false;
        }
        std::string formula;
        if( !translate_attr_string( formula, attr, "formula" )) {
          cerr << "EdmlParser: missing formula for Formula " << name <<  endl;
          return false;
        }
        Formula form;
        form.formula = new TF1(name.c_str(),formula.c_str());
    // define the parameters (should be at least 1)
        for(unsigned ipar=0;ipar<maxnpar;ipar++){
          std::map<std::string, std::string >::const_iterator itr = attr.find( parnames[ipar] );
          if( attr.end() != itr ) {
    // parameter must be defined as a variable
            std::map<std::string,Variable>::const_iterator jpar = _variables.find(itr->second);
            if(jpar == _variables.end() ){
              cerr << "Parameter " << itr->first << " variable " << itr->second << " not found" << endl;
              return false;
            }
            form.parameters.push_back(itr->second);
          }
        }
    // at least 1 variable (parameter) must be defined
        if(form.parameters.size()>0){
          _formulas[name] = form;
        } else {
          cerr << "Formula " << name << " has no variables" << endl;
          return false;
        }
      } else {
        cerr << "EdmlParser:unknown definition " << node << endl;
        return false;
      }
    } while( xmlengine.move2next( noMoreNodes ));
  }
  return noMoreNodes;
}

bool
EdmlParser::translate_loop( EdmlXMLEngine&     xmlengine,
  const std::map<std::string,std::string>& attr,
  Loop& myloop) {
  // get the loop parameters
  if( !translate_loop_attr(attr,myloop)){
    cerr << "Error parsing loop " << endl;
    return false;
  }
//  cout << "found loop of variable " << myloop.loopvar << " with " <<  myloop.loopcount
//    <<" steps starting  with " << myloop.loopstart << " with step " << myloop.loopstep << endl;
  bool noMoreLines;
  if( xmlengine.move2child( noMoreLines )) {
    do {
      const char* ltype = xmlengine.getNodeName( );
      if( !ltype ) return false;
      std::string tstring(ltype);
      std::map<std::string, std::string > lattr;
      if( !xmlengine.getNodeAttr( lattr )) return false;
      myloop.lines.push_back(make_pair(tstring,lattr));
//      cout << "found loop line " << ltype << endl;
  // if this line defines a loop, process it (recursively)
      if(tstring == "loop"){
//        cout << "found nested loop " << endl;
        Loop nestloop;
        bool nloop = translate_loop(xmlengine,lattr,nestloop);
        if(!nloop){
          cerr << "nested loop processing failed " << endl;
          return nloop;
        } else {
  // nest the loops
          myloop.nestedloops.push_back(nestloop);
        }
      }
  //
    } while( xmlengine.move2next( noMoreLines ));
  }
  return true;
}

bool
EdmlParser::execute_loop(const Loop& myloop, Detector& detector,
  Measures& measures,
  const std::string& vname) {
  int loopval = myloop.loopstart;
// find the variable for this loop
  std::map<std::string,Variable>::iterator iloopvar = _variables.find(myloop.loopvar);
  if(iloopvar == _variables.end()){
    cerr << "Loop variable " << myloop.loopvar << " not defined" << endl;
    return false;
  }
  for(int iloop=0;iloop<myloop.loopcount;iloop++){
//    cout << "Executing loop over variable " << myloop.loopvar << " step " << iloop << " with value " << loopval << endl;
    // update the variable for this loop step
    iloopvar->second._ival = loopval;
    // reset nested loop counter
    int inest(0);
    // loop over the lines for this step
    std::vector<std::pair<std::string,std::map<std::string,std::string> > >::const_iterator iline = myloop.lines.begin();
    // legal lines inside a loop are elements or
    while(iline != myloop.lines.end()){
    // test for embedded loops
      if(iline->first == "loop"){
        if(myloop.nestedloops.size() <= inest){
          cerr << "Nested loop mismatch" << endl;
          return false;
        }
      // execute recursively
        execute_loop(myloop.nestedloops[inest++],detector, measures, vname);
      } else {
    // assume everything else is an element
        translate_elem(detector,measures,vname,iline->first.c_str(),iline->second);
      }
      iline++;
    }
    // increment
    loopval += myloop.loopstep;
  }
  return true;
}
