//
//  Model of the target as a diffuse gas
//  Structured text file constructor should have contmin: rmin, rmax, zpos, halfzlength, density (gm/cm^3)
//
#ifndef TrackToy_Detector_Target_hh
#define TrackToy_Detector_Target_hh
#include "TrackToy/Detector/HollowCylinder.hh"
#include "TrackToy/General/TrajUtilities.hh"
#include "TrackToy/Detector/EStar.hh"
#include "TrackToy/General/ELossDistributions.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "TRandom3.h"
#include <string>
#include <vector>
#include <iostream>

namespace TrackToy {
  class Target {
    public:
      using TimeRanges = std::vector<KinKal::TimeRange>;
      enum Material{unknown=-1,Al=1,Ti=2};
      Target(MatEnv::MatDBInfo const& matdbinfo, std::string const& tgtfile,TRandom& tr); // construct from a structured text file
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
      Material matname_;
      const MatEnv::DetMaterial* mat_;
      EStar estar_; // energystar table
      HollowCylinder cyl_;
      double density_;
      double minpath_; // minimum pathlength
      TRandom& tr_; // random number generator
  };

  template<class PKTRAJ> bool Target::extendTrajectory(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj, TimeRanges& intersections,double tol) const {
    using KinKal::VEC3;
    using KinKal::TimeRange;
    static unsigned moyalterms_(20); // number of terms in Moyal expansion
    intersections.clear();
    static double pfactor = 0.001*density_/mat_->density(); // unit conversion cm->mm, and scale for the effective density
    //    std::cout << "density factor" << pfactor << std::endl;
    // extend to the  of the target or exiting the BField (backwards)
    bool retval = extendZ(pktraj,bfield, cyl_.zmax(), tol);
    //    cout << "Z target extend " << ztgt << endl;
    if(retval){ // make sure we didn't exit the BField upstream
      // first find the intersections.
      double tstart = pktraj.range().begin();
      double tstep = 3.0/pktraj.speed(tstart);
      TimeRange trange = cyl_.intersect(pktraj,tstart,tstep);
      while( (!trange.null()) && trange.end() < pktraj.range().end()) {
        double speed = pktraj.speed(trange.mid());
        double plen = std::max(trange.range()*speed, minpath_);
        // require a physical minimum
        intersections.push_back(trange);
        double energy = pktraj.energy(trange.mid());
        auto momvec = pktraj.momentum3(trange.mid());
        double mom = momvec.R();
        // to get physical results, scale the path
//        double estarde = electronEnergyLoss(energy-pktraj.mass(),plen);
        plen *= pfactor;
        double demean = mat_->energyLoss(mom,plen,pktraj.mass());
        double derms = mat_->energyLossRMS(mom,plen,pktraj.mass());
        MoyalDist edist(MoyalDist::MeanRMS(fabs(demean),derms),moyalterms_);
        double ionloss = edist.sample(tr_.Uniform(0.0,1.0));
        // add radiative energy loss.  note we have to convert to cm!!!
        double radFrac = mat_->radiationFraction(trange.range())/10;
        BremssLoss bLoss;
        double bremloss = bLoss.sampleSSPGamma(energy,radFrac);
        // delta energy loss
        DeltaRayLoss dLoss(mat_, energy,plen/10, pktraj.mass());
        double dloss = dLoss.sampleDRL();
        double totloss = ionloss + bremloss + dloss;
        //        std::cout << "Target Ionization eloss = " << ionloss << " Delta eloss " << dloss << " rad eloss "  << bremloss << " tot " << totloss << std::endl;
        energy -= totloss;
//        std::cout << "Target demean " << demean << " derms " << derms << " estar de " << estarde << " totlos " << totloss << std::endl;
        // scattering
        double scatterRMS = mat_->scatterAngleRMS(mom,plen,pktraj.mass());
        //          std::cout << "scatterRMS " << scatterRMS << std::endl;
        // generate random momentum scatter
        VEC3 phidir = VEC3(momvec.Y(),-momvec.X(),0.0).Unit();
        double momtan = tan(momvec.Theta());
        VEC3 thedir = VEC3(momvec.X()/momtan,momvec.Y()/momtan,-momvec.Z()*momtan).Unit();
        auto dmom = tr_.Gaus(0.0,scatterRMS)*mom*phidir + tr_.Gaus(0.0,scatterRMS)*mom*thedir;
        //          std::cout << "phidot " << momvec.Dot(dmom) << " dmag = " << mom - (momvec + dmom).R() << std::endl;
        retval = updateEnergy(pktraj,trange.mid(),energy,dmom);
        if(!retval)break;
        trange = cyl_.intersect(pktraj,trange.end(),tstep);
      }
    }
    return retval;
  }
}
#endif
