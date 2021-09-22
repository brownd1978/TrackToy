//
//  interface to the NIST EStar table, to interpolate between values
//
#ifndef TrackToy_Detector_EStar_hh
#define TrackToy_Detector_EStar_hh
#include <string>
#include <vector>
namespace TrackToy {
  class EStar {
    public:
      // construct from a dump of the table (all 7 columns)
      EStar(const char* tablefile);
      // interpolated values of the tables
      double dEIonization(double energy) const;
      double dERadiation(double energy) const;
      double dETotal(double energy) const;
      double rangeCDSA(double energy) const;
      double radiationYield(double energy) const;
      double densityEffect(double energy) const;
    private:
      std::vector<double> energy_, dEIon_, dERad_, dETot_, range_, radyield_, denseff_; // table
      // find bounding indices; return value says interpolation is needed
      bool findRange(double energy, size_t& imin, size_t& imax) const;
      double interpolate(double energy, std::vector<double> const& table) const;
  };
}
#endif
