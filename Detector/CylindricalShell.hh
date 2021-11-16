//
// Detector element modeled as a 2-dimensional cylindrial shell
// Structured text file constructor should have contain: radius, rhalf, zpos, halfzlength
//
#ifndef TrackToy_Detector_CylindricalShell_hh
#define TrackToy_Detector_CylindricalShell_hh
#include "KinKal/General/TimeRange.hh"
#include "TrackToy/General/TrajUtilities.hh"
#include <string>
#include <vector>
namespace TrackToy {
  using TimeRanges = std::vector<KinKal::TimeRange>;
  class CylindricalShell {
    public:
      CylindricalShell(): radius_(-1.0), rhalf_(-1.0), zpos_(0.0), zhalf_(-1.0) {}
      CylindricalShell(double radius, double rhalf, double zpos, double zhalf) : radius_(radius), rhalf_(rhalf), zpos_(zpos), zhalf_(zhalf) {}
      double radius() const { return radius_;}
      double rhalf() const { return rhalf_;}
      double zmin() const { return zpos_ - zhalf_;}
      double zmax() const { return zpos_ + zhalf_;}
      double zpos() const { return zpos_;}
      double zhalf() const { return zhalf_;}
      // find intersections of a trajectory with this cylinder.  Return the time ranges when the
      // trajectory crosses the shell
      template<class KTRAJ> void intersect(KTRAJ const& ktraj, TimeRanges& tranges,double tstart, double tstep) const;
    private:
      double radius_, rhalf_, zpos_, zhalf_;
  };

  template<class PKTRAJ> void CylindricalShell::intersect(PKTRAJ const& pktraj, TimeRanges& tranges, double tstart, double tstep) const {
    // define boundary times, assuming constant velocity
    tranges.clear();
    double ttest = tstart;
    auto pos = pktraj.position3(ttest);
    double dr = pos.Rho() - radius();
    double olddr = dr;
    while(ttest < pktraj.range().end()){
      //      cout << "particle enters at " << pos << endl;
      auto vel = pktraj.velocity(ttest);
      double dz = fabs(tstep*vel.Z());
      if(pos.Z() > zmin()-dz && pos.Z()< zmax()+dz ){
        // small steps while we're in the z range
        ttest+= tstep;
        auto oldpos = pos;
        pos = pktraj.position3(ttest);
        dr = pos.Rho() - radius();
        if(olddr*dr < 0
            && ( (pos.Z() > zmin() && pos.Z() < zmax()) ||
              (oldpos.Z() > zmin() && oldpos.Z() < zmax()) ) ) {
          // we've crossed the shell.  Interpolate to the exact crossing
          double tx = ttest - tstep*fabs(dr/(dr-olddr));
          // compute the crossing time range
          auto vel = pktraj.velocity(tx);
          double dt = rhalf_/vel.Rho();
          tranges.push_back(KinKal::TimeRange(tx-dt,tx+dt));
        }
        olddr = dr;
      } else {
        // step by piece until we are heading the right direction
        vel = pktraj.velocity(ttest);
        double dt = (zpos() - pos.Z())/vel.Z();
        while(dt < 0.0 && ttest < pktraj.range().end()){
          // advance to the end of this piece
          ttest = pktraj.nearestPiece(ttest).range().end()+tstep;
          pos = pktraj.position3(ttest);
          vel = pktraj.velocity(ttest);
          dt = (zpos() - pos.Z())/vel.Z();
        }
        // now advance to a piece that is within 1 step of the z range
        dt = (vel.Z() > 0) ?  (zmin() - pos.Z())/vel.Z() : (zmax() - pos.Z())/vel.Z();
        while( dt > tstep && ttest < pktraj.range().end()){
          ttest = std::min(ttest+dt, pktraj.nearestPiece(ttest).range().end()+tstep);
          pos = pktraj.position3(ttest);
          vel = pktraj.velocity(ttest);
          dt = (vel.Z() > 0) ?  (zmin() - pos.Z())/vel.Z() : (zmax() - pos.Z())/vel.Z();
        }
      }
    }
  }
}
#endif

