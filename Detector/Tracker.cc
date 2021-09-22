//
//  Implementation of tracker
//
#include "TrackToy/Detector/Tracker.hh"
#include "TrackToy/General/FileFinder.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>
namespace TrackToy {
  Tracker::Tracker(std::string const& tgtfile): ncells_(0), cellArea_(-1.0), cellDensity_(-1.0), activeDensity_(-1.0), gain_(-1.0) {
    FileFinder filefinder;
    std::string fullfile = filefinder.fullFile(tgtfile);
    std::string line;
    static std::string comment("#");
    std::ifstream tgt_stream(fullfile);
    double rmin(-1.0), rmax(-1.0), zpos(-1.0), zhalf(-1.0), density(-1.0);
    unsigned orient(0);
    while (std::getline(tgt_stream, line)) {
      // skip comments and blank lines
      if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
        // strip leading whitespace
        line = std::regex_replace(line, std::regex("^ +"), "");
        std::istringstream iss(line);
        iss >> rmin >> rmax >> zpos >> zhalf >> density >> ncells_ >> orient >> cellArea_ >> activeDensity_ >> gain_;
        break;
      }
    }
    orientation_ = (CellOrientation)orient;
    cyl_ = HollowCylinder(rmin,rmax,zpos,zhalf,density);
    cellDensity_ = (float)ncells_/cyl_.volume();
  }
}


