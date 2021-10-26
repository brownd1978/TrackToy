//
//  Moyal distribution generation
//
#include <cmath>
namespace TrackToy {
// class for Moyal distribution
  class Moyal {
    public:
      Moyal(double mu, double sigma) : mu_(mu), sigma_(sigma) {}
      double mu() const { return mu_; }
      double sigma() const { return sigma_; }
      double mean() const;
      double variance() const;
      double RMS() const;
      // generate a moyal distribution.  Input is assumed random over the interval [0-1];
      double sample(double randval) const;
      double PDF(double x) const;
    private:
      double mu_;
      double sigma_;
  };
}
