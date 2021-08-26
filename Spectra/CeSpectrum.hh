//
//  Model of the Ce spectrum including radiative effects, see https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.051301
//
namespace TrackToy {
  struct CeSpectrumParams {
    double eMax_; // endpoint energy
    CeSpectrumParams(double emax) : eMax_(emax) {}
  };

  class CeSpectrum {
    public:
      CeSpectrum(CeSpectrumParams const& params);
      // relative probability to observe this energy (per MeV)
      double rate(double energy) const;
      // integral over a range
      double integral(double elow, double ehigh, unsigned nstep=1000 ) const;
    private:
      CeSpectrumParams params_;
      double norm_; // normalization to rate/MeV
  };


}
