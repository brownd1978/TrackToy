//
//  Simple class to find a file in TrackToy
//
#ifndef TrackToy_General_FileFinder_hh
#define TrackToy_General_FileFinder_hh
#include <string>
#include <cstdlib>
#include <stdexcept>
namespace TrackToy {
  class FileFinder {
    public:
      FileFinder(std::string envname="TRACKTOY_SOURCE_DIR") : envname_(envname) {
        const char* env =  std::getenv(envname_.c_str());
        if(!env) throw std::invalid_argument("TRACKTOY_SOURCE_DIR not set: did you forget to source setup.sh?");
        sourcedir_ = std::string(env);
      }
      std::string fullFile(std::string const& file) const {
        return sourcedir_+"/"+file;
      }
    private:
      std::string envname_, sourcedir_;
  };
}
#endif
