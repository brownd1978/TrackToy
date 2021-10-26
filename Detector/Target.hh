//
//  Model of the target as a diffuse gas
//  Structured text file constructor should have contmin: rmin, rmax, zpos, halfzlength, density (gm/cm^3)
//
#ifndef TrackToy_Detector_Target_hh
#define TrackToy_Detector_Target_hh
#include "TrackToy/Detector/HollowCylinder.hh"
#include "TrackToy/General/TrajUtilities.hh"
#include "TrackToy/Detector/EStar.hh"
#include "KinKal/General/TimeRange.hh"
#include <string>
#include <vector>

namespace TrackToy {
  class Target {
    public:
      enum Material{unknown=-1,Al=1,Ti=2};
      Target(std::string const& tgtfile); // construct from a structured text file
      double density() const { return density_;} // gm/cm^3
      auto cylinder() const { return cyl_; }
      // find the  energy loss (mean and RMS) from estar for a path through the target
      double electronEnergyLoss(double ke, double pathlen) const;
      double protonEnergyLoss(double ke, double pathlen) const; // TODO
      std::string material() const;
      // extrapolate the trajectory through the target and update it for the energy loss on exit
      template<class PKTRAJ> double updateTrajectory(PKTRAJ& pktraj) const;
    private:
      Material mat_;
      EStar estar_; // energystar table
      HollowCylinder cyl_;
      double density_;
  };

  template<class PKTRAJ> double Target::updateTrajectory(PKTRAJ& pktraj) const {
    double retval(0.0);
    std::vector<KinKal::TimeRange> tranges;
    static double tstep(0.01);
    cyl_.intersect(pktraj,tranges,tstep);
    if(tranges.size() > 0){
      double energy = pktraj.energy(tranges.front().begin());
      double speed = pktraj.speed(tranges.front().begin());
      for (auto const& range : tranges) {
        double pathlen = range.range()*speed;
        // should check for particle type FIXME!
        double de = electronEnergyLoss(energy-pktraj.mass(),pathlen);
        energy += de;
        retval += de;
        if(energy < pktraj.mass())break;
      }
      if(energy > pktraj.mass()){
 //       double endtime = tranges.back().end();
//        updateEnergy(pktraj,endtime,retval);
      }
    }
    return retval;
  }
}
#endif
