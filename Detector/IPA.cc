//
//  Implementation of IPA
//
#include "KinKal/Detector/MaterialXing.hh"
#include "TrackToy/Detector/IPA.hh"
#include "TrackToy/General/FileFinder.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>
#include <string>
namespace TrackToy {
  IPA::IPA(MatEnv::MatDBInfo const& matdbinfo,std::string const& tgtfile) : type_(unknown) {
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
        // first get type and material
        if(type_ == unknown){
          int type;
          std::string material;
          iss >> type >> material;
          type_ = (IPAType)type;
          // lookup material
          mat_ = matdbinfo.findDetMaterial(material);
          if(mat_ == 0){
            std::string errmsg = std::string("Invalid Material ")+material;
            throw std::invalid_argument(errmsg);
          }
          // then geometry
        } else if (type_ == cylinder) {
          std::cout << "line=" << line << std::endl;
          double radius, rhalf, zpos, zhalf;
          iss >> radius >> rhalf >> zpos >> zhalf;
          if(radius < 0.0 || rhalf < 0.0 || zhalf < 0.0)throw std::invalid_argument("Invalid CylindricalShell parameters");\
            cyl_ = CylindricalShell(radius,rhalf,zpos,zhalf);
        } else if (type_ == propeller) {
          throw std::invalid_argument("Propeller IPA not currently supported");\
        }
      }
    }
    switch (type_ ) {
      default:
      case IPA::cylinder:
        std::cout << "Cylindrical IPA with radius " << cyl_.radius() << " thickness " << 2*cyl_.rhalf() << " Zmid " << cyl_.zpos() << " Length " << 2*cyl_.zhalf();
        break;
      case IPA::propeller:
        std::cout << "Propeller IPA";
        break;
    }
    std::cout << " of material " << mat_->name() << " constructed from file " << fullfile << std::endl;
  }
}

