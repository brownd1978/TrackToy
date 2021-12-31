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
  CeMinusSpectrum::CeMinusSpectrum(CeMinusSpectrumParams const& params) : Spectrum(Spectrum::CeMinus, params.EEnd_, params.EEnd_),params_(params) {
  }

  double CeMinusSpectrum::rate(double energy) const {
    if(energy == params_.EEnd_)
      return 1.0;
    else
      return 0.0;
  }

  double CeMinusSpectrum::integral(double elow, double ehigh, unsigned nsteps)const {
    if(elow < params_.EEnd_ && ehigh > params_.EEnd_)
      return 1.0;
    else
      return 0.0;
  }

  double CeMinusSpectrum::sample(double prob) const {
    return params_.EEnd_;
  }
}
