//
//  Model of the CeMinus spectrum including radiative effects, see https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.051301
//
#ifndef TrackToy_Spectra_CeMinusSpectrum_hh
#define TrackToy_Spectra_CeMinusSpectrum_hh
namespace TrackToy {
  struct CeMinusSpectrumParams {
    double eMax_; // endpoint energy
    CeMinusSpectrumParams(double emax) : eMax_(emax) {}
  };

  class CeMinusSpectrum {
    public:
      CeMinusSpectrum(CeMinusSpectrumParams const& params);
      // relative probability to observe this energy (per MeV)
      double rate(double energy) const;
      // integral over a range
      double integral(double elow, double ehigh, unsigned nstep=10000 ) const;
      auto const& params() const { return params_; }
    private:
      CeMinusSpectrumParams params_;
      double norm_; // normalization to rate/MeV
  };
}
#endif
