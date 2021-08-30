//
//  Model of the CeMinus spectrum including radiative effects, see https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.051301
//
#include "TrackToy/Spectra/CeMinusSpectrum.hh"
#include "TrackToy/General/PhysicalConstants.h"
#include <cmath>
#include <iostream>

using CLHEP::fine_structure_const;
using CLHEP::electron_mass_c2;

namespace TrackToy {
  CeMinusSpectrum::CeMinusSpectrum(CeMinusSpectrumParams const& params) : params_(params), norm_(1.0) {
    // normalize
    norm_ = 1.0/integral(params_.eMax_,0.0);
  }

  double CeMinusSpectrum::rate(double energy) const {
    double retval(0.0);
    if(energy > 0.0 && energy < params_.eMax_)
      retval = norm_*(fine_structure_const/(2*M_PI))*(log(4*energy*energy/electron_mass_c2/electron_mass_c2)-2.)*((energy*energy+params_.eMax_*params_.eMax_)/(params_.eMax_-energy));
    return retval;
  }

  double CeMinusSpectrum::integral(double elow, double ehigh, unsigned nsteps)const {
    double estep = (ehigh-elow)/nsteps;
    double emin = std::max(0.0,elow);
    double emax = std::min(ehigh,params_.eMax_);
    double sum = 0.0;
    for( unsigned istep=0; istep<nsteps;++istep){
      double energy = emin+(istep+0.5)*estep;
      sum += rate(energy);
    }
    return estep*sum/(emax-emin);
  }

}
