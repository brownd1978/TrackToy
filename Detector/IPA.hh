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
        double energy = pktraj.energy(intersections.front().begin());
        for (auto const& ipainter : intersections) {
          double mom = pktraj.momentum(ipainter.mid());
          double plen = pktraj.speed(ipainter.mid())*ipainter.range();
          // Moyal dist. models ionization loss
          double demean = mat_->energyLoss(mom,plen,pktraj.mass());
          double derms = mat_->energyLossRMS(mom,plen,pktraj.mass());
          // model ionization energy loss using a Moyal distribution
          MoyalDist edist(MoyalDist::MeanRMS(fabs(demean), derms),10);
          double ionloss = edist.sample(tr_.Uniform(0.0,1.0));
          // add radiative energy loss.  note we have to convert to cm!!!
          double radFrac = mat_->radiationFraction(ipainter.range()/10.0);
          KinKal::BremssLoss bLoss;
          double bremloss = bLoss.sampleSSPGamma(energy,radFrac);
          // delta energy loss
          KinKal::DeltaRayLoss dLoss(mat_, mom,plen, pktraj.mass());
          double dloss = dLoss.sampleDRL();
          double totloss = ionloss + bremloss + dloss;
//          std::cout << "IPA Ionization eloss = " << ionloss << " Delta eloss " << dloss << " rad eloss "  << bremloss << " tot " << totloss << std::endl;
          double oldenergy = energy;
          energy -= totloss;
          std::cout << "old energy " << oldenergy << " new energy " << energy << std::endl;
          retval = updateEnergy(pktraj,ipainter.end(),energy);
          if(!retval)break;
        }
      }
    }
    return retval;
  }
}
#endif

