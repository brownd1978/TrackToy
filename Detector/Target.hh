//
//  Model of the target as a diffuse gas
//  Structured text file constructor should have contmin: rmin, rmax, zpos, halfzlength, density (gm/cm^3)
//
#ifndef TrackToy_Detector_Target_hh
#define TrackToy_Detector_Target_hh
#include "TrackToy/Detector/HollowCylinder.hh"
#include "TrackToy/General/TrajUtilities.hh"
#include "TrackToy/Detector/EStar.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/General/TimeRange.hh"
#include <string>
#include <vector>
#include <iostream>

namespace TrackToy {
  class Target {
    public:
      using TimeRanges = std::vector<KinKal::TimeRange>;
      enum Material{unknown=-1,Al=1,Ti=2};
      Target(std::string const& tgtfile); // construct from a structured text file
      double density() const { return density_;} // gm/cm^3
      auto cylinder() const { return cyl_; }
      // find the  energy loss (mean and RMS) from estar for a path through the target
      double electronEnergyLoss(double ke, double pathlen) const;
      double protonEnergyLoss(double ke, double pathlen) const; // not yet implemented TODO
      std::string material() const;
      // extend a trajectory through the target.  Return value specifies if the particle continues downsream (true) or stops in the target or exits the field upstream (false)
      template<class PKTRAJ> bool extendTrajectory(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj,TimeRanges& intersections,double tol=1e-4) const;
      void print(std::ostream& os) const;
    private:
      Material mat_;
      EStar estar_; // energystar table
      HollowCylinder cyl_;
      double density_;
  };

  template<class PKTRAJ> bool Target::extendTrajectory(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj, TimeRanges& intersections,double tol) const {
    bool retval(false);
    intersections.clear();
    // extend to the  of the target or exiting the BField (backwards)
    retval = extendZ(pktraj,bfield, cyl_.zmax(), tol);
    //    cout << "Z target extend " << ztgt << endl;
    if(retval){ // make sure we didn't exit the BField upstream
      // first find the intersections.
      double tstart = pktraj.range().begin();
      double tstep = 3.0/pktraj.speed(tstart);
      cyl_.intersect(pktraj,intersections,tstart,tstep);
      if(intersections.size() > 0){
        double energy = pktraj.energy(intersections.front().begin());
        double speed = pktraj.speed(intersections.front().begin());
        for (auto const& range : intersections) {
          double pathlen = range.range()*speed;
          // should check for particle type FIXME!
          double de = electronEnergyLoss(energy-pktraj.mass(),pathlen);
          energy -= de;
        }
        retval = updateEnergy(pktraj,intersections.back().end(),energy);
      }
    }
    return retval;
  }
}
#endif
