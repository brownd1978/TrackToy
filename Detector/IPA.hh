//
//  Inner proton absorber.  Currently a thin cylinder
//
#ifndef TrackToy_Detector_IPA_hh
#define TrackToy_Detector_IPA_hh
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/General/TimeRange.hh"
#include "TrackToy/General/Moyal.hh"
#include "TrackToy/Detector/CylindricalShell.hh"
namespace TrackToy {
  class IPA {
    public:
      enum IPAType { unknown=-1, cylinder=1, propeller };
      IPA(MatEnv::MatDBInfo const& matdbinfo,std::string const& file);
      auto const& cyl() const { return cyl_; }
      auto const& material() const { return *mat_; }
      auto type() const { return type_; }
      // find the Moyal distribution describing the eloss function for a given intersection of a trajectory with this cylinder.
      template<class PKTRAJ> Moyal energyLoss(PKTRAJ const& pktraj, KinKal::TimeRange const& trange) const;
      template<class PKTRAJ> bool updateTrajectory(PKTRAJ& pktraj, TimeRanges& intersections) const;
    private:
      IPAType type_;
      CylindricalShell cyl_;
      const MatEnv::DetMaterial* mat_;

  };

  template<class PKTRAJ> bool IPA::updateTrajectory(PKTRAJ& pktraj, TimeRanges& intersections) const {
    bool retval(false);
    // first find the intersections.
    double tstart = pktraj.range().begin();
    static double tstep(0.01);
    cyl_.intersect(pktraj,intersections,tstart,tstep);
    if(intersections.size() > 0){
      double energy = pktraj.energy(intersections.front().begin());
      double speed = pktraj.speed(intersections.front().begin());
      for (auto const& range : intersections) {
	double pathlen = range.range()*speed;
	// should check for particle type FIXME!
	double de = electronEnergyLoss(energy-pktraj.mass(),pathlen);
	energy += de;
      }
      retval = updateEnergy(pktraj,intersections.back().end(),energy);
    }
    return retval;
  }

  template<class PKTRAJ> Moyal IPA::energyLoss(PKTRAJ const& pktraj, KinKal::TimeRange const& trange) const {
    // sample the trajectory at the middle of the range
    auto const& ktraj = pktraj.nearestPiece(trange.mid());
    double beta = ktraj.beta();
    double mom = ktraj.momentum();
    double plen = ktraj.speed(trange.mid())*trange.range();
    double xi = mat_->eloss_xi(beta,plen);
    double dp = mat_->energyLossMPV(mom,plen,ktraj.mass());
    return Moyal(dp,xi);
  }
}
#endif

