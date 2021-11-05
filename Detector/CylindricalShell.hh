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
    // find the front
    double ttest = ztime(pktraj,tstart,zmin());
    auto pos = pktraj.position3(ttest);
    double dr = pos.Rho() - radius();
    double olddr = dr;
    // add test for looping TODO
    while(ttest < pktraj.range().end() && pos.Z() < zmax()){
      // step
      ttest += tstep;
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
    }
  }
}
#endif

