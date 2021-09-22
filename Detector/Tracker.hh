//
// Cylindrical tracker
//
#ifndef TrackToy_Detector_Tracker_hh
#define TrackToy_Detector_Tracker_hh
#include "TrackToy/Detector/HollowCylinder.hh"
#include <string>
#include <vector>
namespace TrackToy {
  class Tracker {
    public:
      enum CellOrientation{radial=0,axial};
      Tracker(std::string const& tfile);
      auto const& cylinder() const { return cyl_; }
      double cellArea() const { return cellArea_; }
      unsigned nCells() const { return ncells_; }
      // compute # of cells for a given trajectory
      template<class KTRAJ> double nCells(KTRAJ const& ktraj, TimeRanges const& tranges)  const;
    private:
      HollowCylinder cyl_; // geometric form of the tracker
      CellOrientation orientation_; // orientation of the cells
      unsigned ncells_; // number of cells
      double cellArea_; // average x-sectional area of the cells (not including track projection)
      double cellDensity_; // cells per mm^3
      double activeDensity_;  // gm/cm^3 of active material
      double gain_; // gain factor (C/MeV)
  };

  template<class KTRAJ> double Tracker::nCells(KTRAJ const& ktraj, TimeRanges const& tranges) const {
    double retval(0.0);
    if(tranges.size()>0){
      // assume constant speed
      double speed = ktraj.speed(tranges.front().begin());
      double path(0.0);
      for(auto const& range : tranges) path += range.range()*speed;
      retval = cellDensity_*path*cellArea_;
    }
    return retval;
  }
}
#endif

