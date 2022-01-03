//
//  Model of the CeMinus spectrum including radiative effects, see https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.051301
//
#ifndef TrackToy_Spectra_CeMinusSpectrum_hh
#define TrackToy_Spectra_CeMinusSpectrum_hh
#include "TrackToy/Spectra/Spectrum.hh"
#include <limits>
 namespace TrackToy {
  struct CeMinusSpectrumParams {
    double EEnd_; // endpoint energy
    CeMinusSpectrumParams(double emax) : EEnd_(emax) {}
  };

  class CeMinusSpectrum : public Spectrum {
    public:
      CeMinusSpectrum(CeMinusSpectrumParams const& params);
      double rate(double energy) const override;
      double integral(double elow, double ehigh, unsigned nstep=100 ) const override;
      double sample(double prob) const override;
      auto const& params() const { return params_; }
    private:
      CeMinusSpectrumParams params_;
  };
}
#endif
