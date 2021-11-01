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
    // compute the time tolerance based on the speed.  Require 1mm precision (good enough for IPA
    double ttol = 1.0/pktraj.speed(pktraj.range().begin());

//
//      if(targetinters.size() > 0) tstart = targetinters.back().end();
//      // extend through ipa
//      extendZ(pktraj,axfield, axfield.zMin(), ipa.cylinder().zmax(), tol);
//      // find intersections with ipa
//      ipa.cylinder().intersect(pktraj,ipainters,tstart,tstep);
//      nipa_ = ipainters.size();
//      nipa->Fill(nipa_);
//      ipade_ = 0.0;
//      for(auto const& ipainter : ipainters) {
//        auto eloss = ipa.energyLoss(pktraj,ipainter);
//        double de = eloss.mean();
//        ipade->Fill(-de);
//        double der = eloss.sample(tr_.Uniform(0.0,1.0)); // currently broken, FIXME
//        ipader->Fill(-der);
//        //        ipade_ += der;
//        ipade_ += de;
//      }
//      ipades->Fill(-ipade_);
//      bool ipacont(true);
//      if(ipainters.size() > 0){
//        double energy = pktraj.energy(tstart) + ipade_;
//        ipacont = updateEnergy(pktraj,ipainters.back().end(),energy);
//        tstart = ipainters.back().end();
//      }

// record the end of the previous extension; this is where new extensions start
    double tstart = pktraj.back().range().begin();
    double energy = pktraj.energy(tstart);
    // extend through the IPA or exiting the BField (backwards)
    double ztgt = extendZ(pktraj,bfield, bfield.zMin(), cyl_.zmax(), ttol);
    //    cout << "Z target extend " << ztgt << endl;
    if(ztgt > cyl_.zmin()){ // make sure we didn't exit the BField upstream
      // first find the intersections.
      static double tstep(0.01);
      cyl_.intersect(pktraj,intersections,tstart,tstep);
      if(intersections.size() > 0){
        for (auto const& ipainter : intersections) {
          auto eloss = energyLoss(pktraj,ipainter);
          double de = eloss.mean(); // should sample the distribution FIXME
//        double de = eloss.sample(tr_.Uniform(0.0,1.0)); // currently broken, FIXME
          energy -= de;
        }
        retval = updateEnergy(pktraj,intersections.back().end(),energy);
      }
    }
    return retval;
  }
}
#endif

