//
//  Implementation of HollowCylinder
//
#include "TrackToy/Detector/HollowCylinder.hh"
#include "TrackToy/General/FileFinder.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>
namespace TrackToy {
  HollowCylinder::HollowCylinder(std::string const& tgtfile): HollowCylinder() {
    FileFinder filefinder;
    std::string fullfile = filefinder.fullFile(tgtfile);
    std::string line;
    static std::string comment("#");
    std::ifstream tgt_stream(fullfile);
    while (std::getline(tgt_stream, line)) {
      // skip comments and blank lines
      if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
        // strip leading whitespace
        line = std::regex_replace(line, std::regex("^ +"), "");
        std::istringstream iss(line);
        iss >> rmin_ >> rmax_ >> zpos_ >> zhalf_ >> density_;
      }
    }
    if(rmin_ < 0.0 || rmax_ < 0.0 || zhalf_ < 0.0 || density_ < 0.0)throw std::invalid_argument("Invalid HollowCylinder parameters");\
      std::cout << "Read HollowCylinder parameters from file " << fullfile << std::endl;
  }

}
