//
//  Inner proton absorber.  Currently a thin cylinder
//
#ifndef TrackToy_Detector_IPA_hh
#define TrackToy_Detector_IPA_hh
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/General/TimeRange.hh"
#include "TrackToy/General/Moyal.hh"
#include "TrackToy/General/TrajUtilities.hh"
#include "TrackToy/Detector/CylindricalShell.hh"
namespace TrackToy {
  class IPA {
    public:
      enum IPAType { unknown=-1, cylindrical=1, propeller };
      IPA(MatEnv::MatDBInfo const& matdbinfo,std::string const& file);
      auto const& cylinder() const { return cyl_; }
      auto const& material() const { return *mat_; }
      auto type() const { return type_; }
      // find the Moyal distribution describing the eloss function for a given intersection of a trajectory with this cylinder.
      template<class PKTRAJ> Moyal energyLoss(PKTRAJ const& pktraj, KinKal::TimeRange const& trange) const;
      // extend a trajectory through the IPA.  Return value specifies if the particle continues downsream (true) or stops in the target or exits the field upstream (false)
      template<class PKTRAJ> bool extendTrajectory(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj,TimeRanges& intersections) const;
      void print(std::ostream& os) const;
    private:
      IPAType type_;
      CylindricalShell cyl_;
      const MatEnv::DetMaterial* mat_;
  };

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

  template<class PKTRAJ> bool IPA::extendTrajectory(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj, TimeRanges& intersections) const {
    bool retval(false);
    intersections.clear();
    // compute the time tolerance based on the speed.
    double ttol = 3.0/pktraj.speed(pktraj.range().begin());
    // record the end of the previous extension; this is where new extensions start
    double tstart = pktraj.back().range().begin();
    double energy = pktraj.energy(tstart);
    // extend through the IPA or exiting the BField (backwards)
    retval = extendZ(pktraj,bfield, cyl_.zmax(), ttol);
//   std::cout << "IPA extend " << retval << std::endl;
    if(retval){
//      auto pstart = pktraj.position3(tstart);
//      std::cout << "zstart " << pstart.Z() << " Z extend " << zipa << std::endl;
      // first find the intersections.
      static double tstep(0.01);
      cyl_.intersect(pktraj,intersections,tstart,tstep);
      if(intersections.size() > 0){
        for (auto const& ipainter : intersections) {
          auto eloss = energyLoss(pktraj,ipainter);
          double de = eloss.mean(); // should sample the distribution FIXME
//        double de = eloss.sample(tr_.Uniform(0.0,1.0)); // currently broken, FIXME
          energy += de;
        }
        retval = updateEnergy(pktraj,intersections.back().end(),energy);
//        std::cout << "IPA update " << retval << std::endl;
      }
    }
    return retval;
  }
}
#endif

