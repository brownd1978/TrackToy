//
//  Implementation of calorimeter
//
#include "TrackToy/Detector/Calorimeter.hh"
#include "TrackToy/General/FileFinder.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>

namespace TrackToy {
  Calorimeter::Calorimeter(std::string const& calofile): vprop_(-1.0), tres_(-1.0), pres_(-1.0), minpath_(-1.0)
  {
    FileFinder filefinder;
    std::string fullfile = filefinder.fullFile(calofile);
    std::string line;
    static std::string comment("#");
    std::ifstream tgt_stream(fullfile,std::ios_base::in);
    if(tgt_stream.fail()){
      std::string errmsg = std::string("File doesn't exist" )+ fullfile;
      throw std::invalid_argument(errmsg.c_str());
    }
    while (std::getline(tgt_stream, line)) {
      // skip comments and blank lines
      if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
        // strip leading whitespace
        line = std::regex_replace(line, std::regex("^ +"), "");
        std::istringstream iss(line);
        // first, disk 1 geometry
        double rmin, rmax, zpos, zhalf;
        if(disks_[0].rmin() < 0.0){
          iss >> rmin >> rmax >> zpos >> zhalf;
          disks_[0] = HollowCylinder(rmin,rmax,zpos,zhalf);
          // disk 2 geometry
        } else if (disks_[1].rmin() < 0.0){
          iss >> rmin >> rmax >> zpos >> zhalf;
          disks_[2] = HollowCylinder(rmin,rmax,zpos,zhalf);
        } else {
          iss >> vprop_ >> tres_ >> pres_ >> minpath_;
        }
      }
    }
  }

  void Calorimeter::print(std::ostream& os ) const {
    std::cout << "Calorimeter  disk 1 " << disk(0) << " disk 2 " << disk(1) << std::endl;
    std::cout << " Vprop " << vProp() << " time resolution " << timeResolution()
    << " position resolution " << positionResolution() << " min hit path " << minPath() << std::endl;
  }
}


