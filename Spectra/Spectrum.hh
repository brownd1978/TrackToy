//
//  interface for kinetic energy spectra
//
#ifndef TrackToy_Spectra_Spectrum_hh
#define TrackToy_Spectra_Spectrum_hh
#include <string>
#include <vector>
#include <limits>
namespace TrackToy {
  class Spectrum {
    public:
      enum SpectrumType {DIO=0, CeMinus, CeEndpoint};
      Spectrum(SpectrumType type, double emin=-1., double emax = std::numeric_limits<float>::max()) : type_(type), eMin_(emin), eMax_(emax) {}
      // rate to observe this kinetic energy (per MeV), relative to the total integral
      virtual double rate(double energy) const =0;
      // sample the distribution; the input variable must be in the interval [0,1]
      virtual double sample(double prob) const =0;
      // integral rate to observe kinetic energies within the range
      virtual double integral(double elow, double ehigh, unsigned nstep=1000 ) const=0;
      // basic infrastructure
      virtual double eMin() const { return eMin_; }
      virtual double eMax() const { return eMax_; }
      virtual double normalization() const { return norm_; }
      SpectrumType type() const { return type_; }
    protected:
      SpectrumType type_; // type of spectrum
      double eMin_, eMax_; // lower and upper energies of this spectrum
      double norm_; // normalization to the integral rate for this interval
  };
}
#endif

