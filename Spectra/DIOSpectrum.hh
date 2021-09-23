//
//  Model of the captured mu- Decay In Orbit (DIO) spectrum including radiative effects, see https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.051301
//
#ifndef TrackToy_Spectra_DIOSpectrum_hh
#define TrackToy_Spectra_DIOSpectrum_hh
#include <string>
#include <vector>
#include <limits>
namespace TrackToy {
  class DIOSpectrum {
    public:
      DIOSpectrum(const char* spectrumfile,double emin=-1.0,double emax = std::numeric_limits<float>::max());
      // relative probability to observe this energy (per MeV)
      double rate(double energy) const;
      // integrate over a range
      double integral(double elow, double ehigh, unsigned nstep=1000 ) const;
      // sample the distribution; the input variable must be in the interval [0,1]
      double sample(double prob) const;
      double eMin() const { return eMin_; }
      double eMax() const { return eMax_; }
      double normalization() const { return norm_; }
    private:
      std::vector<double> rate_; // tabulated rates
      std::vector<double> cdf_; // CDF (for random sampling) with the same energy binning
      double eMin_, eMax_, eStep_; // lower and upper energies of the rate table
      double norm_; // normalization to rate/stopped muon/MeV
  };
}
#endif
