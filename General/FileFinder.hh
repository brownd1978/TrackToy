//
//  Simple class to find a file in TrackToy
//
#ifndef TrackToy_General_FileFinder_hh
#define TrackToy_General_FileFinder_hh
#include <string>
#include <cstdlib>
namespace TrackToy {
  class FileFinder {
    public:
      FileFinder(std::string envname="TRACKTOY_SOURCE_DIR") : envname_(envname) {
        sourcedir_ = std::getenv(envname_.c_str());
      }
      std::string fullFile(std::string const& file) const {
        return sourcedir_+"/"+file;
      }
    private:
      std::string envname_, sourcedir_;
  };
}
#endif
