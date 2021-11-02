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
  Tracker::Tracker(MatEnv::MatDBInfo const& matdbinfo,std::string const& tgtfile):
    ncells_(0), cellDensity_(-1.0), smat_(0)
  {
    FileFinder filefinder;
    std::string fullfile = filefinder.fullFile(tgtfile);
    std::string line;
    static std::string comment("#");
    std::ifstream tgt_stream(fullfile);
    double rmin(-1.0), rmax(-1.0), zpos(-1.0), zhalf(-1.0);
    double rcell(-1.0), lcell(-1.0), rwire(-1.0), wthick(-1.0);
    unsigned orient(0);
    while (std::getline(tgt_stream, line)) {
      // skip comments and blank lines
      if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
        // strip leading whitespace
        line = std::regex_replace(line, std::regex("^ +"), "");
        std::istringstream iss(line);
        // first, overall geometry
        if(rmin < 0.0){
          iss >> rmin >> rmax >> zpos >> zhalf;
          cyl_ = HollowCylinder(rmin,rmax,zpos,zhalf);
        } else if(rcell < 0.0) {
          // then cell description
          iss >> orient >> ncells_ >> rcell >> lcell >> wthick >> rwire ;
          std::cout << "ncells " << ncells_  << " rcell " << rcell << " lcell " << lcell << " wthick " << wthick << " rwire " << rwire << std::endl;
          if(rcell<0.0)throw std::invalid_argument("Invalid cell parameters");
          smat_ = new KinKal::StrawMaterial(matdbinfo, rcell, wthick, rwire);
          orientation_ = (CellOrientation)orient;
          // cell density (transverse to the cell)
          cellDensity_ = ncells_*2.0*rcell*lcell/cyl_.volume();
          // compute the total mass of 1 cell
          double cmass = smat_->wallMaterial().density()*2.0*M_PI*rcell*lcell*wthick
            + smat_->gasMaterial().density()*M_PI*rcell*rcell*lcell
            + smat_->wireMaterial().density()*M_PI*rwire*rwire*lcell;
          density_ = 1000.0*cmass*ncells_/cyl_.volume(); // convert to gm/cm^3
        }
      }
    }
  }

  void Tracker::print(std::ostream& os ) const {
    std::cout << "Tracker between " << cyl_.zmin() << " and " << cyl_.zmax() << " rmin " << cyl_.rmin() << " rmax " << cyl_.rmax()
    << " material density (gm/cm^3)" << density_ << " cell density (cells/mm) " << cellDensity_ << " with " << ncells_ << " cells oriented ";
    if(orientation_ == azimuthal)
      std::cout << " azimuthally " << std::endl;
    else
      std::cout << " axially " << std::endl;
  }
}


