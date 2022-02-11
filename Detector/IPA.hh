//
//  Inner proton absorber.  Currently a thin cylinder
//
#ifndef TrackToy_Detector_IPA_hh
#define TrackToy_Detector_IPA_hh
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/MatEnv/ELossDistributions.hh"
#include "KinKal/General/TimeRange.hh"
#include "TrackToy/General/TrajUtilities.hh"
#include "TrackToy/Detector/CylindricalShell.hh"
#include "TRandom3.h"
namespace TrackToy {
  class IPA {
    public:
      enum IPAType { unknown=-1, cylindrical=1, propeller };
      IPA(MatEnv::MatDBInfo const& matdbinfo,std::string const& file,TRandom& tr);
      auto const& cylinder() const { return cyl_; }
      auto const& material() const { return *mat_; }
      auto type() const { return type_; }
      // extend a trajectory through the IPA.  Return value specifies if the particle continues downsream (true) or stops in the target or exits the field upstream (false)
      template<class PKTRAJ> bool extendTrajectory(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj,TimeRanges& intersections,double tol=1e-4) const;
      void print(std::ostream& os) const;
    private:
      IPAType type_;
      CylindricalShell cyl_;
      const MatEnv::DetMaterial* mat_;
      TRandom& tr_; // random number generator
  };


  template<class PKTRAJ> bool IPA::extendTrajectory(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj, TimeRanges& intersections,double tol) const {
    using KinKal::MoyalDist;
    bool retval(false);
    intersections.clear();
    // record the end of the previous extension; this is where new extensions start
    double tstart = pktraj.back().range().begin();
    double energy = pktraj.energy(tstart);
    // extend through the IPA or exiting the BField (backwards)
    retval = extendZ(pktraj,bfield, cyl_.zmax(), tol);
//   std::cout << "IPA extend " << retval << std::endl;
    if(retval){
//      auto pstart = pktraj.position3(tstart);
//      std::cout << "zstart " << pstart.Z() << " Z extend " << zipa << std::endl;
      // first find the intersections.
      static double tstep(0.01);
      cyl_.intersect(pktraj,intersections,tstart,tstep);
      if(intersections.size() > 0){
        for (auto const& ipainter : intersections) {
          double mom = pktraj.momentum(ipainter.mid());
          double plen = pktraj.speed(ipainter.mid())*ipainter.range();
          // Moyal dist. models ionization loss
          double demean = mat_->energyLoss(mom,plen,pktraj.mass());
          double derms = mat_->energyLossRMS(mom,plen,pktraj.mass());
          MoyalDist edist(MoyalDist::MeanRMS(demean, derms),10);
          double de = edist.sample(tr_.Uniform(0.0,1.0));
          // add radiative energy loss: model as an expontential with the same mean as the ionization (ie critical energy).
          // this is particle/energy specific
//          if(pktraj.beta(ipainter.mid()) > 0.99)de -= tr_.Exp(fabs(demean));
          //          std::cout << "IPA de " << de << std::endl;
          energy += de;
        }
        retval = updateEnergy(pktraj,intersections.back().end(),energy);
      }
    }
    return retval;
  }
}
#endif

