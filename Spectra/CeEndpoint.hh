//
// trivial spectrum of CeEndpoint
//
#ifndef TrackToy_Spectra_CeEndpoint_hh
#define TrackToy_Spectra_CeEndpoint_hh
#include "TrackToy/Spectra/Spectrum.hh"
#include <limits>
 namespace TrackToy {

  class CeEndpoint : public Spectrum {
    public:
      CeEndpoint(double EEnd) : Spectrum(Spectrum::CeEndpoint, EEnd-1.0, EEnd+1.0), EEnd_(EEnd) { norm_ = 1.0; }
      double rate(double energy) const override { if(energy == EEnd_) return 1.0; else return 0.0; }
      double integral(double elow, double ehigh, unsigned nstep=10000 ) const override { if(EEnd_ > elow && EEnd_ < ehigh) return 1.0; else return 0.0; }
      double sample(double prob) const override { return EEnd_; }
      double endpointEnergy() const { return EEnd_; }
    private:
      double EEnd_; // endpoint energy
  };
}
#endif
