#include "TrackToy/Detector/Target.hh"
#include "TrackToy/General/FileFinder.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>

namespace TrackToy {

  Target::Target(std::string const& tgtfile): mat_(unknown), density_(-1.0) {
    FileFinder filefinder;
    std::string fullfile = filefinder.fullFile(tgtfile);
    std::string line;
    static std::string comment("#");
    std::ifstream tgt_stream(fullfile);

    std::string efile_al("Data/EStar_Al.dat");

    double rmin, rmax, zpos, zhalf;
    std::string material;
    while (std::getline(tgt_stream, line)) {
      // skip comments and blank lines
      if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
        // strip leading whitespace
        line = std::regex_replace(line, std::regex("^ +"), "");
        std::istringstream iss(line);
        if(mat_ == unknown){
          iss >> material >> density_;
          if(material == "Al"){
            mat_ = Al;
            // create the EStar from this
            estar_ = EStar(efile_al);
          } else {
            std::string errmsg = std::string("Invalid Material ")+material;
            throw std::invalid_argument(errmsg);
          }
        } else {
          // get cylinder paramters
          iss >> rmin >> rmax >> zpos >> zhalf;
        }
      }
    }
    if(rmin < 0.0 || rmax < 0.0 || zhalf < 0.0 || density_ < 0.0)throw std::invalid_argument("Invalid Target parameters");
    cyl_ = HollowCylinder(rmin,rmax,zpos,zhalf);
    std::cout << "Read Target from file " << fullfile << std::endl;
  }

  std::string Target::material() const {
    std::string retval("unknown");
    switch (mat_) {
    case Al:
      retval = "Aluminum";
      break;
    default:
      break;
    }
    return retval;
  }

  double Target::electronEnergyLoss(double ke, double pathlen) const {
    return estar_.dEIonization(ke)*density()*pathlen/10.0; // only ionization energy loss is relevant for thin material
    // note this corrects for density being in gm/cm^3 (estar table is relative to that)
  }

}
