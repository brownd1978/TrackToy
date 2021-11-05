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
    //      cout << "particle enters at " << pos << endl;
    if(pos.Z() < zmin() || pos.Z() > zmax()) {
      // advance to where we enter the volume going the right way
      auto tmin = ztime(pktraj,tstart,zmin());
      auto tmax = ztime(pktraj,tstart,zmax());
      ttest = std::min(tmin,tmax)+1.0e-6; //  step a tiny amount
      pos = pktraj.position3(ttest);
    }
    // if we start inside, setup the first range
    bool inside = isInside(pos);
    bool oldinside = inside;
    bool crosses(false);
    if(inside) tranges.push_back(KinKal::TimeRange(ttest,ttest));
    while(ttest < pktraj.range().end() && pos.Z() < zmax()){
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
      // if we stepped outside the z range, see if the particle loops back
      if(pos.Z() < zmin() ) {
        // advance to where we enter the volume going the right way
        auto tmin = ztime(pktraj,ttest,zmin())-1.0e-6;
        ttest = std::max(ttest,tmin);
        pos = pktraj.position3(ttest);
        inside = isInside(pos);
        oldinside = inside;
      }
    }
    // finish the last range
    if(inside)tranges.back().end() = std::min(ttest,pktraj.range().end());
  }
}
#endif
