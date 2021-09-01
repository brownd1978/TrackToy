//
//  Compute the range and other properties of a muon of a givem momentum in a given material description file from the PDG
//  https://pdg.lbl.gov/2021/AtomicNuclearProperties/
//  all non-tabular names must be commented (#)
//
#ifndef TrackToy_General_MuonRange_hh
#define TrackToy_General_MuonRange_hh
#include <vector>
namespace TrackToy {
// Range is reported as CSDA g/cm^2
  struct RangeData {
    RangeData(double energy, double momentum, double range) : energy_(energy), momentum_(momentum), range_(range) {}
    double energy_, momentum_, range_;
  };

  class MuonRange {
    public:
      MuonRange(const char* datafile, double density); // density in gm/cm^3
// find the range in mm given muon kinetic energy in MeV
      double rangeEnergy(double energy) const;
// find the range in mm given muon momentum in MeV/c
      double rangeMomentum(double mom) const;
    private:
      double density_; // gm/cm^3
      std::vector<RangeData> range_;
  };
}
#endif
