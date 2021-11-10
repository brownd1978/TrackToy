//
// Detector element with hollow right cylinder geometry
//
#ifndef TrackToy_Detector_HollowCylinder_hh
#define TrackToy_Detector_HollowCylinder_hh
#include <string>
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/Vectors.hh"
#include "TrackToy/General/TrajUtilities.hh"
#include <vector>
#include <algorithm>
namespace TrackToy {
  using TimeRanges = std::vector<KinKal::TimeRange>;
  class HollowCylinder {
    public:
      HollowCylinder(): rmin_(-1.0), rmax_(-1.0), zpos_(0.0), zhalf_(-1.0) {}
      HollowCylinder(double rmin, double rmax, double zpos, double zhalf) : rmin_(rmin), rmax_(rmax), zpos_(zpos), zhalf_(zhalf) {}
      double rmin() const { return rmin_;}
      double rmax() const { return rmax_;}
      double zmin() const { return zpos_ - zhalf_;}
      double zmax() const { return zpos_ + zhalf_;}
      double zpos() const { return zpos_;}
      double zhalf() const { return zhalf_;}
      double length() const { return 2.0*zhalf_; }
      double area() const { return M_PI*(rmax_*rmax_ - rmin_*rmin_); }
      double volume() const { return length()*area(); } // mm^3
      // find intersections of a trajectory with this cylinder.  Return the time ranges in which the
      // trajectory is inside the physical volume
      template<class PKTRAJ> void intersect(PKTRAJ const& pktraj, TimeRanges& tranges, double tstart, double tstep) const;
    private:
      double rmin_, rmax_, zpos_, zhalf_;
      bool isInside(KinKal::VEC3 const& pos) const {
        double rho = pos.Rho();
        return rho > rmin() && rho < rmax() && pos.Z() > zmin() && pos.Z() < zmax();
      }
  };

  template<class PKTRAJ> void HollowCylinder::intersect(PKTRAJ const& pktraj, TimeRanges& tranges, double tstart, double tstep) const {
    tranges.clear();
    // see if we are starting inside
    double ttest = tstart;
    auto pos = pktraj.position3(ttest);
    bool inside = isInside(pos);
    if(inside) tranges.push_back(KinKal::TimeRange(ttest,ttest));
    while(ttest < pktraj.range().end()){
      //      cout << "particle enters at " << pos << endl;
      auto vel = pktraj.velocity(ttest);
      double dz = fabs(tstep*vel.Z());
      if(pos.Z() > zmin()-dz && pos.Z()< zmax()+dz ){
        // small steps while we're in the z range
        bool oldinside = inside;
        bool crosses(false);
        ttest+= tstep;
        pos = pktraj.position3(ttest);
        oldinside = inside;
        inside = isInside(pos);
        crosses = oldinside != inside;
        //          cout << "ttest " << ttest << " pos " << pos << endl;
        if(crosses){
          if(oldinside){
            // finish this range
            tranges.back().end() = ttest;
          }
          else {
            // entering: create the range
            tranges.push_back(KinKal::TimeRange(ttest,ttest));
          }
        }
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
        inside = isInside(pos);
        if(inside) tranges.push_back(KinKal::TimeRange(ttest,ttest));
      }
    }
    // finish the last range
    if(inside)tranges.back().end() = std::min(ttest,pktraj.range().end());
  }
}
#endif
