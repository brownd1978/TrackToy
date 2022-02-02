//
// twin disk calorimeter
//
#ifndef TrackToy_Detector_Calorimeter_hh
#define TrackToy_Detector_Calorimeter_hh
#include "TrackToy/Detector/HollowCylinder.hh"
#include "TrackToy/General/TrajUtilities.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Trajectory/Line.hh"
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
      Calorimeter(std::string const& tfile,TRandom& tr);
      auto const& disk(size_t idisk) const { return disks_[idisk]; }
      auto const& vProp() const { return vprop_; }
      double showerMax() const { return shmax_; }
      double timeResolution() const { return tres_; }
      double positionResolution() const { return pres_; }
      double minPath() const { return minpath_; }
      void print(std::ostream& os) const;
      template<class KTRAJ> unsigned simulateHits(KinKal::BFieldMap const& bfield,
          KinKal::ParticleTrajectory<KTRAJ>& mctraj,
          std::vector<std::shared_ptr<KinKal::Hit<KTRAJ>>>& hits, double tol) const;
      template <class KTRAJ> void simulateHit(KinKal::ParticleTrajectory<KTRAJ> const& mctraj,
          KinKal::TimeRange const& trange, std::vector<std::shared_ptr<KinKal::Hit<KTRAJ>>>& hits) const;
    private:
      std::array<HollowCylinder,2> disks_; // 2 calorimeter disks
      KinKal::VEC3 vprop_; // light propagation velocity
      double tres_, pres_; // time and (transverse) position resolution
      double shmax_; // shower max depth
      double minpath_; // minimum path to generate a hit
      TRandom& tr_; // random number generator
  };
  template<class KTRAJ> unsigned Calorimeter::simulateHits(KinKal::BFieldMap const& bfield,
      KinKal::ParticleTrajectory<KTRAJ>& mctraj,
      std::vector<std::shared_ptr<KinKal::Hit<KTRAJ>>>& hits, double tol) const {
    unsigned retval(0);
   // extend through the first disk
    extendZ(mctraj,bfield,disk(0).zmax(),tol);
    double tstart = ztime(mctraj,mctraj.back().range().begin(),disk(0).zmin()-10.0);
    double speed = mctraj.speed(tstart);
    double tstep = 0.1*pres_/speed; // set step to fraction of transverse size
    // find intersections with the first disk
    std::vector<KinKal::TimeRange> tinters;
    disk(0).intersect(mctraj,tinters,tstart,tstep);
    if(tinters.size()== 0) {
    // go to the second disk
      extendZ(mctraj,bfield,disk(0).zmax(),tol);
      tstart = ztime(mctraj,mctraj.back().range().begin(),disk(1).zmin()-10.0);
      disk(1).intersect(mctraj,tinters,tstart,tstep);
    }
    if(tinters.size() > 0 && tinters.front().range()*speed > minpath_){
      simulateHit(mctraj,tinters.front(),hits);
      retval ++;
    }
    return retval;
  }

  template <class KTRAJ> void Calorimeter::simulateHit(KinKal::ParticleTrajectory<KTRAJ> const& mctraj,
    KinKal::TimeRange const& trange, std::vector<std::shared_ptr<KinKal::Hit<KTRAJ>>>& hits) const {
    using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,KinKal::Line>;
    using SCINTHIT = KinKal::ScintHit<KTRAJ>;
    using KinKal::Line;
    using KinKal::CAHint;
    // put the hit at showermax, or the end of the range (with a buffer)
    double speed = mctraj.speed(trange.begin());
    double tshmax = std::min(trange.begin() + shmax_/speed, trange.end()-0.1);
    auto shmaxpos = mctraj.position3(tshmax);
    auto hitpos = shmaxpos;
    // figure out which disk
    size_t idisk = (hitpos.Z() > disks_[0].zmax()) ? 1 : 0;
    // smear the transverse position
    hitpos.SetX(tr_.Gaus(hitpos.X(),pres_));
    hitpos.SetY(tr_.Gaus(hitpos.Y(),pres_));
    // set the z position to the sensor plane (back of the disk)
    hitpos.SetZ(disks_[idisk].zmax());
    double dz = disks_[idisk].zmax()-shmaxpos.Z();
    // set the measurement time to correspond to the light propagation from showermax_, smeared by the resolution
    double tmeas = tr_.Gaus(tshmax+dz/vprop_.Z(),tres_);
    Line hline(hitpos,tmeas,vprop_,disks_[idisk].length());
    // create the hit and add it; the hit has no material, since the original particle doesn't exist after showering so we don't need to simulate its interactions
    CAHint tphint(tmeas,tmeas);
    static double tprec(1e-8); // TPOCA precision
    PTCA tp(mctraj,hline,tphint,tprec);
    hits.push_back(std::make_shared<SCINTHIT>(tp, tres_*tres_, pres_*pres_,tprec));
  }
}

#endif

