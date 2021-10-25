#include "TrackToy/General/Moyal.hh"
namespace TrackToy {

  // inverse erf isn't defined in std; use an approximation here
  float erfinvf (float a);

  double Moyal::mean() const {
    constexpr static double moyalmeanfactor = 0.57721566490153286 + M_LN2 ; //approximate Euler-Mascheroni (also known as gamma) constant (0.577...), see https://mathworld.wolfram.com/Euler-MascheroniConstant.html, added to log(2). This sum is used for the calculation of the closed-form Moyal mean below
    return mu_ + sigma_ * moyalmeanfactor; //formula from https://reference.wolfram.com/language/ref/MoyalDistribution.html, see end of file for more information
  }

  double Moyal::variance() const {
    static double const pisq2(0.5*M_PI*M_PI);
    return pisq2*sigma_*sigma_;
  }
  double Moyal::RMS() const {
    static double const pi_over_sqrt2(M_PI/M_SQRT2);
    return pi_over_sqrt2*sigma_;
  }

  // generate a moyal distribution.  Input is assumed random over the interval [0-1];
  double Moyal::sample(double randval) const {
    float erfval = 1.0-2*randval;
    float inverfc = erfinvf(erfval);
    static double const sqrt2(sqrt(2.0));
    double moyalval = mu_ - 2.0*sigma_*log(sqrt2*inverfc);
    return moyalval;
  }

  double Moyal::PDF(double x) const {
    static double const invsqrt2pi(1.0/sqrt(2.0*M_PI));
    double diff = (x-mu_)/sigma_;
    return invsqrt2pi*exp( -0.5*(diff + exp(diff)));
  }

  // see https://stackoverflow.com/questions/27229371/inverse-error-function-in-c/49743348#49743348
 float erfinvf (float a)
  {
    float p, r, t;
    t = std::fmaf (a, 0.0f - a, 1.0f);
    t = logf(t);
    if (fabsf(t) > 6.125f) { // maximum ulp error = 2.35793
      p =              3.03697567e-10f; //  0x1.4deb44p-32
      p = std::fmaf (p, t,  2.93243101e-8f); //  0x1.f7c9aep-26
      p = std::fmaf (p, t,  1.22150334e-6f); //  0x1.47e512p-20
      p = std::fmaf (p, t,  2.84108955e-5f); //  0x1.dca7dep-16
      p = std::fmaf (p, t,  3.93552968e-4f); //  0x1.9cab92p-12
      p = std::fmaf (p, t,  3.02698812e-3f); //  0x1.8cc0dep-9
      p = std::fmaf (p, t,  4.83185798e-3f); //  0x1.3ca920p-8
      p = std::fmaf (p, t, -2.64646143e-1f); // -0x1.0eff66p-2
      p = std::fmaf (p, t,  8.40016484e-1f); //  0x1.ae16a4p-1
    } else { // maximum ulp error = 2.35002
      p =              5.43877832e-9f;  //  0x1.75c000p-28
      p = std::fmaf (p, t,  1.43285448e-7f); //  0x1.33b402p-23
      p = std::fmaf (p, t,  1.22774793e-6f); //  0x1.499232p-20
      p = std::fmaf (p, t,  1.12963626e-7f); //  0x1.e52cd2p-24
      p = std::fmaf (p, t, -5.61530760e-5f); // -0x1.d70bd0p-15
      p = std::fmaf (p, t, -1.47697632e-4f); // -0x1.35be90p-13
      p = std::fmaf (p, t,  2.31468678e-3f); //  0x1.2f6400p-9
      p = std::fmaf (p, t,  1.15392581e-2f); //  0x1.7a1e50p-7
      p = std::fmaf (p, t, -2.32015476e-1f); // -0x1.db2aeep-3
      p = std::fmaf (p, t,  8.86226892e-1f); //  0x1.c5bf88p-1
    }
    r = a * p;
    return r;
  }
}
