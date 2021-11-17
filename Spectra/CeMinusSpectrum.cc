//
//  Model of the CeMinus spectrum including radiative effects, see https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.051301
//
#include "TrackToy/Spectra/CeMinusSpectrum.hh"
#include "KinKal/General/PhysicalConstants.h"
#include <cmath>
#include <iostream>

using CLHEP::fine_structure_const;
using CLHEP::electron_mass_c2;

namespace TrackToy {
  CeMinusSpectrum::CeMinusSpectrum(CeMinusSpectrumParams const& params, double emin, double emax) : Spectrum(Spectrum::CeMinus, std::max(0.0,emin),std::min(params.EEnd_,emax)),
  params_(params) {
    // normalize
    norm_ = 1.0/integral(eMin(),eMax());
    unsigned nsteps(100000);
    eStep_ = (eMax_ - eMin_)/(nsteps-1);
    // create the CDF
    cdf_.reserve(nsteps);
    for(size_t istep=0;istep < nsteps; ++istep){
      double pcdf = cdf_.size() > 0 ? cdf_.back() : 0.0;
      double energy = eMin_ + (istep)*eStep_;
      double rv = rate(energy)*norm_;
      cdf_.push_back(pcdf + rv*eStep_);
    }
  }

  double CeMinusSpectrum::rate(double energy) const {
    double retval(0.0);
    if(energy > 0.0 && energy < params_.EEnd_)
      retval = (fine_structure_const/(2*M_PI))*(log(4*energy*energy/electron_mass_c2/electron_mass_c2)-2.)*((energy*energy+params_.EEnd_*params_.EEnd_)/(params_.EEnd_-energy));
    return retval;
  }

  double CeMinusSpectrum::integral(double elow, double ehigh, unsigned nsteps)const {
    double estep = (ehigh-elow)/nsteps;
    double emin = std::max(eMin(),elow);
    double emax = std::min(eMax(),params_.EEnd_);
    double sum = 0.0;
    for( unsigned istep=0; istep<nsteps;++istep){
      double energy = emin+(istep+0.5)*estep;
      sum += rate(energy);
    }
    return estep*sum/(emax-emin);
  }

  double CeMinusSpectrum::sample(double prob) const {
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


}
