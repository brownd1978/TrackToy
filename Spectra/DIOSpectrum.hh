//
//  Model of the captured mu- Decay In Orbit (DIO) spectrum including radiative effects, see https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.051301
//
#ifndef TrackToy_Spectra_DIOSpectrum_hh
#define TrackToy_Spectra_DIOSpectrum_hh
#include <string>
#include <vector>
namespace TrackToy {
  class DIOSpectrum {
    public:
      DIOSpectrum(const char* spectrumfile);
      // relative probability to observe this energy (per MeV)
      double rate(double energy) const;
      // integral over a range
      double integral(double elow, double ehigh, unsigned nstep=1000 ) const;
    private:
      std::vector<double> rate_; // tabulated rates
      double eMin_, eMax_, eStep_; // lower and upper energies of the rate table
      double norm_; // normalization to rate/stopped muon/MeV
  };
}
#endif
