//
//  Model of the captured mu- Decay In Orbit (DIO) spectrum including radiative effects, see https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.051301
//
#include "TrackToy/Spectra/DIOSpectrum.hh"
#include "TrackToy/General/FileFinder.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <iterator>
#include <cmath>


namespace TrackToy {
  DIOSpectrum::DIOSpectrum(const char* spectrumfile, double emin, double emax) : Spectrum(Spectrum::DIO, emin,emax)  {
    FileFinder filefinder;
    std::string fullfile = filefinder.fullFile(spectrumfile);
    // read the file
    std::ifstream spectrum_stream(fullfile);
    std::string line;
    bool first(true);
    static std::string comment("#");
    double energy, oldenergy(-1.0), val;
    while (std::getline(spectrum_stream, line)) {
      // skip comments and blank lines
      if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
        std::istringstream iss(line);
        while ((iss >> energy >> val )) {
          if(energy > emin && energy < emax) {
            if(first){
              first = false;
              eMin_ = energy;
            } else {
              // check spectrum is uniform
              double de = energy-oldenergy;
              size_t test = rint((energy-eMin_)/de);
              if(test != rate_.size()) std::cout << "Error: spetrum isn't uniform!" <<
                " energy " << energy << " de " << de << " rate " <<  val << " expected index " << test << " size " << rate_.size() << std::endl;
            }
            rate_.push_back(val);
            oldenergy = energy;
          }
          eMax_ = oldenergy;
        }
      }
    }
    eStep_ = (eMax_ - eMin_)/(rate_.size()-1);
    // create the CDF
    norm_ = 1.0/integral(eMin_, eMax_,rate_.size());
    cdf_.reserve(rate_.size());
    for(size_t irate=0;irate < rate_.size(); ++irate){
      double pcdf = cdf_.size() > 0 ? cdf_.back() : 0.0;
 //      double energy = eMin_ + 0.5*eStep_;
      double energy = eMin_ + (irate+0.5)*eStep_;
      double rv = rate(energy)*norm_;
      cdf_.push_back(pcdf + rv*eStep_);
    }
    std::cout << "CDF diff " << cdf_.back()-1.0 << std::endl;
//    if(cdf_.back()-1.0>1e-5)throw std::range_error("Invalid CDF");
  }

  double DIOSpectrum::rate(double energy) const {
    double retval(0.0);
    if(energy < eMin_){
      // linear extrapolation; assume 0 rate at 0 energy (?)
      if(eMin_ > 0.0){
        double drdE = rate_[0]/(eMin_);
        retval = std::max(0.0,rate_[0]+drdE*(energy-eMin_));
      } else
        retval = 0.0;
    } else if(energy > eMax_){
      // linear extrapolation
      double drdE = (rate_[rate_.size()-1] - rate_[rate_.size()-2])/eStep_;
      retval = std::max(0.0,rate_[rate_.size()-1]+drdE*(energy-eMax_));
    } else {
      // linear interpolation
      size_t imin = static_cast<size_t>( std::floor((energy-eMin_)/eStep_) );
      double drdE = (rate_[imin+1] - rate_[imin])/eStep_;
      retval = rate_[imin] + drdE*(energy - eMin_ - eStep_*imin);
    }
    return retval;
  }

  double DIOSpectrum::sample(double prob) const {
    if(prob < 0.0 || prob > 1.0)throw std::invalid_argument("Invalid probability");
    double retval(eMax_);
    // find the nearest bin
    size_t ibin(0);
    while(ibin < cdf_.size() && prob > cdf_[ibin])++ibin;
    if(ibin > 0 && ibin < cdf_.size()-1){
      size_t jbin = ibin-1;
      retval = eMin_ + (jbin + (prob-cdf_[jbin])/(cdf_[ibin]-cdf_[jbin]))*eStep_;
    } else if (ibin == 0) {
      retval = eMin_;
    } else
      retval = eMax_;
    return retval;
  }

  double DIOSpectrum::integral(double elow, double ehigh, unsigned nsteps)const {
    double emin = std::max(0.0,eMin_);
    double emax = std::min(ehigh,eMax_);
    double estep = (emax-emin)/nsteps;
    double sum = 0.0;
    for( unsigned istep=0; istep<nsteps;++istep){
      double energy = emin+(istep+0.5)*estep;
      sum += rate(energy);
    }
    return estep*sum;
  }
}
