/* Class: EdmlParser
 *
 * A class for parsing for EDML TrackToy configurations.
 */
#ifndef EdmlParser_HH
#define EdmlParser_HH

#include "TrackToy/Detector/EdmlXMLEngine.hh"
#include <iostream>
#include <string>
#include <vector>
#include <map>

class EdmlXMLEngine;
template <class a > class EdmlAnyTypeDict;

struct Formula{
  std::vector<std::string> parameters;
  TF1* formula;
};

struct Variable{
  enum vtype{intvar=0,floatvar,stringvar};
  vtype _type;
  int _ival;
  float _fval;
  std::string _sval;
};

///  The parser class for EDML configurations
/**
 * The class is a stateful fron-end interface to the actual parser.
 * It will follows parsing requests and cache results in its data
 * members.
 */
class EdmlParser {

  private:

    /* These methods and are prohibited because of the instance-level cache
     * which can't be trivially cloned in its current implementation.
     */
    EdmlParser( const EdmlParser& );
    EdmlParser& operator=( const EdmlParser& );

  public:

    /// The c-tor

    EdmlParser( const bool verbose=false );

    /// The d-tor
    /**
     * NOTE: It will destroy all subdetectors (detector elements)
     *       owned by an object.
     */
    virtual ~EdmlParser( );

    /// Parse the specified file
    /**
     * The method will attempt to locate and parse the file and if
     * successfull return 'true'.
     *
     * OWNERSHIP NOTE: the ownership is NOT returned by the method
     */
    bool parse( const char* filename );

    /// Access the names of files parsed
    /**
     * The method will return a list of files used for parsing between
     * the calls to the reset methods.
     */
    void filenames( std::vector<std::string >& theNames ) const;

    /// Access the latest detector description
    /**
     * The method will return 0 if no parsing has been done or
     * if the last parsing attempt failed.
     *
     * OWNERSHIP NOTE: the ownership is NOT returned by the method
     */
    const EdmlDetector* detector( ) const { return _detector; }

    /// Get a list of known measuring devices
    /**
     * Device identifiers (unique names) serve as the dictionary keys.
     * And the values are the corresponding device objects.
     * A detailed information on each class can be found in a documentation
     * section for the below mentioned method:
     *
     * @see method EdmlMeasuringDevice::type()
     */
    void devices( std::vector<EdmlMeasuringDevice >& theDevices ) const;

    /// Get a dictionary of known configuration parameters
    /**
     * Parameter names are keys to the dictionary. The corresponding values are
     * represented pairs of the entry's type and the entry's value exactly as it
     * was found in the input file.
     */
    void config( std::map<std::string,
        std::pair<std::string,
        std::string > >& theParameters ) const;

    /// Dump results of the last successfull (if any) parsing
    /**
     * The output will be printed onto the specified stream and 'true'
     * returned. If the last parsing was either non-successfull, or if
     * it never happend then nothing will be printed onto the specified
     * output stream and the 'false' result will be returned instead.
     */
    bool dump( std::ostream& str = std::cout ) const;

    /// Reset a context of the last operations
    /**
     * A description of the method to be added later.
     */
    void reset( );

  private:
    // translation functions

    bool translate_attr_float( float& value,
        const std::map<std::string, std::string >& attr,
        const char*                                name );
    bool translate_attr_float( EdmlAnyTypeDict<std::string>&              dict,
        const std::map<std::string, std::string >& attr,
        const char*                                name );
    bool translate_attr_float_if_avail( float&                                     value,
        const std::map<std::string, std::string >& attr,
        const char*                                name );
    bool translate_attr_int( int&                                     value,
        const std::map<std::string, std::string >& attr,
        const char*                                name );
    bool translate_attr_int_if_avail( int&                                     value,
        const std::map<std::string, std::string >& attr,
        const char*                                name );
    bool translate_attr_float_if_avail( EdmlAnyTypeDict<std::string>&              dict,
        const std::map<std::string, std::string >& attr,
        const char*                                name );
    bool translate_attr_string( std::string&                               value,
        const std::map<std::string, std::string >& attr,
        const char*                                name );
    bool translate_attr_string( EdmlAnyTypeDict<std::string>&              dict,
        const std::map<std::string, std::string >& attr,
        const char*                                name );
    bool translate_attr_string_if_avail( std::string&                               value,
        const std::map<std::string, std::string >& attr,
        const char*                                name );
    bool translate_attr_string_if_avail( EdmlAnyTypeDict<std::string>&              dict,
        const std::map<std::string, std::string >& attr,
        const char*                                name );

        // TrackToy specifics
    bool translate_elem( Detector&        detector,
        Measures&        measures,
        const std::string& vname, const char* type,const std::map<std::string,std::string>& attr);
    bool translate_volume( Detector&        detector,
        Measures&        measures,
        EdmlXMLEngine&     xmlengine,
        const std::string& vname );
    bool translate_detector( Detector&    detector,
        Measures&    measures,
        EdmlXMLEngine& xmlengine );
    bool translate_config( Config&          config,
        EdmlXMLEngine&     xmlengine,
        const std::string& basepath = "" );
    bool translate( Detector& detector,
        Measures& measures,
        Config&   config,
        const char* filename );

    // translate attribute values using formulas, variables, definitions, or literals
    std::string translate_attr_value( const std::string& value) const;
    // deal with loops
    bool translate_loop_attr(const std::map<std::string,std::string>& attr,Loop& loop);
    bool translate_loop( EdmlXMLEngine&     xmlengine,
        const std::map<std::string,std::string>& attr,Loop& loop);
    bool execute_loop(const Loop& myloop,   Detector&        detector,
        Measures&        measures,
        const std::string& vname);

    bool translate_definitions(  EdmlXMLEngine& xmlengine );

    // Parameters controlling the parser's behavior
    //
    bool _verbose;

    // Cached parameters & results of the last parsing
    //
    std::vector<std::string > _filenames;


    std::map<std::string,std::string> _constants;
    std::map<std::string,Variable> _variables;
    std::map<std::string,Formula> _formulas;

    std::map<std::string,                       // path
      std::pair<std::string,             // type
      std::string > > _config; // value
};

#endif // EdmlParser_HH
