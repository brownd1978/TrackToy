//
//  EStart interface implementation
//
#include "TrackToy/Detector/EStar.hh"
#include "TrackToy/General/FileFinder.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

namespace TrackToy {
  EStar::EStar(std::string const& tablefile) {
    FileFinder filefinder;
    std::string fullfile = filefinder.fullFile(tablefile);
    // read the file
    std::ifstream spectrum_stream(fullfile);
    std::string line;
    static std::string comment("#");
    double energy, dEIon, dERad, dETot, range, radyield, denseff;
    while (std::getline(spectrum_stream, line)) {
      // skip comments and blank lines
      if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
        std::istringstream iss(line);
        while ((iss >> energy >> dEIon >> dERad >> dETot >> range >> radyield >> denseff )) {
          energy_.push_back(energy);
          dEIon_.push_back(dEIon);
          dERad_.push_back(dERad);
          dETot_.push_back(dETot);
          range_.push_back(range);
          radyield_.push_back(radyield);
          denseff_.push_back(denseff);
        }
      }
    }
    std::cout << "Read " << energy_.size() << " rows from EStar table " << fullfile << std::endl;
  }
  double EStar::dEIonization(double energy) const { return interpolate(energy,dEIon_); }
  double EStar::dERadiation(double energy) const{ return interpolate(energy,dERad_); }
  double EStar::dETotal(double energy) const { return interpolate(energy,dETot_); }
  double EStar::rangeCDSA(double energy) const { return interpolate(energy,range_); }
  double EStar::radiationYield(double energy) const { return interpolate(energy,radyield_); }
  double EStar::densityEffect(double energy) const { return interpolate(energy,denseff_); }

  double EStar::interpolate(double energy, std::vector<double>const& table) const {
    size_t imin(0), imax(0);
    if(findRange(energy,imin,imax)){
      double vmin = table[imin];
      double dv = table[imax] - vmin;
      double de = energy_[imax]-energy_[imin];
      return vmin + (energy-energy_[imin])*dv/de;
    } else
      return table[imin];
  }

  bool EStar::findRange(double energy, size_t& imin, size_t& imax) const {
    if(energy <= energy_.front()){
      imin = imax = 0;
      return false;
    } else if(energy >= energy_.back()){
      imin = imax = energy_.size()-1;
      return false;
    } else {
      imax = 0;
      while(energy_[imax] < energy)++imax;
      imin = imax-1;
      return true;
    }
  }
}
