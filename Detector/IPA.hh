//
//  Inner proton absorber.  Currently a thin cylinder
//
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
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
  double elossMPV = mat_->energyLoss(mom,plen,ktraj.mass()); // most probable energy loss = Moyal mu
  return Moyal(elossMPV, mat_->eloss_xi(beta,plen));

  }
}


