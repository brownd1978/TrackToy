//
//  Implementation of CylinderElem
//
#include "TrackToy/Detector/CylinderElem.hh"
#include "TrackToy/General/FileFinder.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>
namespace TrackToy {
  CylinderElem::CylinderElem(std::string const& tgtfile): rmin_(-1.0), rmax_(-1.0), zpos_(0.0), zhalf_(-1.0), density_(-1.0), adensity_(0.0), ncells_(0) {
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
        iss >> rmin_ >> rmax_ >> zpos_ >> zhalf_ >> density_ >> adensity_ >> ncells_;
      }
    }
    if(rmin_ < 0.0 || rmax_ < 0.0 || zhalf_ < 0.0 || density_ < 0.0)throw std::invalid_argument("Invalid CylinderElem parameters");\
      std::cout << "Read CylinderElem parameters from file " << fullfile << std::endl;
  }

}
