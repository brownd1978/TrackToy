//
// twin disk calorimeter
//
#ifndef TrackToy_Detector_Calorimeter_hh
#define TrackToy_Detector_Calorimeter_hh
#include "TrackToy/Detector/HollowCylinder.hh"
#include "TRandom3.h"
#include "Math/VectorUtil.h"
#include "KinKal/Examples/ScintHit.hh"
#include <array>
#include <string>
#include <vector>
#include <iostream>
namespace TrackToy {
  class Calorimeter {
    public:
      Calorimeter(std::string const& tfile);
      auto const& disk(size_t idisk) const { return disks_[idisk]; }
      double vProp() const { return vprop_; }
      double timeResolution() const { return tres_; }
      double positionResolution() const { return pres_; }
      double minPath() const { return minpath_; }
      void print(std::ostream& os) const;
      template<class KTRAJ> void simulateHits(KinKal::BFieldMap const& bfield,
          KinKal::ParticleTrajectory<KTRAJ>& mctraj,
          std::vector<std::shared_ptr<KinKal::Hit<KTRAJ>>>& hits,
          std::vector<KinKal::TimeRange>& tinters, double tol) const;
    private:
      std::array<HollowCylinder,2> disks_; // 2 calorimeter disks
      double vprop_; // light propagation velocity
      double tres_, pres_; // time and (transverse) position resolution
      double minpath_; // minimum path to generate a hit
  };
}

#endif

