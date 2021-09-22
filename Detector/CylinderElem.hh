//
// Detector element with hollow right cylinder geometry
// Structured text file constructor should have contmin: rmin, rmax, zpos, halfzlength, density (gm/cm^3)
//
#ifndef TrackToy_Detector_CylinderElem_hh
#define TrackToy_Detector_CylinderElem_hh
#include <string>
#include "KinKal/General/TimeRange.hh"
#include <vector>
namespace TrackToy {
  using TimeRanges = std::vector<KinKal::TimeRange>;
  class CylinderElem {
    public:
      CylinderElem(std::string const& tgtfile); // construct from a structured text file
      double rmin() const { return rmin_;}
      double rmax() const { return rmax_;}
      double zmin() const { return zpos_ - zhalf_;}
      double zmax() const { return zpos_ + zhalf_;}
      double zpos() const { return zpos_;}
      double zhalf() const { return zhalf_;}
      double density() const { return density_;}
      double activeDensity() const { return adensity_;}
      // find intersections of a trajectory with this cylinder.  Return the time ranges in which the
      // trajectory is inside the physical volume
      template<class KTRAJ> void intersect(KTRAJ const& ktraj, TimeRanges& tranges, double tstep) const;
    private:
      double rmin_, rmax_, zpos_, zhalf_, density_, adensity_;
      unsigned ncells_;
  };

  template<class KTRAJ> void CylinderElem::intersect(KTRAJ const& ktraj, TimeRanges& tranges, double tstep) const {
    // define boundary times, assuming constant velocity
    double tmin = ktraj.zTime(zmin());
    tranges.clear();
    // search for an intersection; start with the entrance
    auto pos = ktraj.position3(tmin);
    double rho = pos.Rho();
    bool inside = rho > rmin() && rho < rmax();
    bool oldinside;
    bool crosses(false);
    //      cout << "particle enters at " << pos << endl;
    double ttest = tmin;
    // if we start inside, setup the first range
    if(inside) tranges.push_back(KinKal::TimeRange(tmin,tmin));
    while(pos.Z() < zmax()){
      ttest += tstep;
      pos = ktraj.position3(ttest);
      double rho = pos.Rho();
      oldinside = inside;
      inside = rho > rmin() && rho < rmax();
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
    if(inside)tranges.back().end() = ttest;
  }
}
#endif

