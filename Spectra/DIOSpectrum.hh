//
//  Model of the captured mu- Decay In Orbit (DIO) spectrum including radiative effects, see https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.051301
//
#ifndef TrackToy_Spectra_DIOSpectrum_hh
#define TrackToy_Spectra_DIOSpectrum_hh
#include "TrackToy/Spectra/Spectrum.hh"
#include <string>
#include <vector>
#include <limits>
namespace TrackToy {
  class DIOSpectrum : public Spectrum {
    public:
      DIOSpectrum(const char* spectrumfile,double emin=-1.0,double emax = std::numeric_limits<float>::max());
      double rate(double energy) const override;
      double integral(double elow, double ehigh, unsigned nstep=1000 ) const override;
      double sample(double prob) const override;
      double eStep() const { return eStep_; }
    private:
      std::vector<double> rate_; // tabulated rates
      std::vector<double> cdf_; // CDF (for random sampling) with the same energy binning
      double eStep_; // energy step
  };
}
#endif
