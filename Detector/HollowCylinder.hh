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
      double volume() const { return 2.0*M_PI*zhalf_*(rmax_*rmax_ - rmin_*rmin_); } // mm^3
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
    bool oldinside = inside;
    bool crosses(false);
    //      cout << "particle enters at " << pos << endl;
    // if we start inside, setup the first range
    if(inside) tranges.push_back(KinKal::TimeRange(ttest,ttest));
    while(ttest < pktraj.range().end()){
    // step
      ttest = timeStep(pktraj,zmin(),zmax(),ttest,tstep);
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
    }
    // finish the last range
    if(inside)tranges.back().end() = std::min(ttest,pktraj.range().end());
  }
}
#endif
