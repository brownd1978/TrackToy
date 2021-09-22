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
      double nCells(double speed, TimeRanges const& tranges) const;
    private:
      HollowCylinder cyl_; // geometric form of the tracker
      CellOrientation orientation_; // orientation of the cells
      unsigned ncells_; // number of cells
      double cellArea_; // average x-sectional area of the cells (not including track projection)
      double cellDensity_; // cells per mm^3
      double activeDensity_;  // gm/cm^3 of active material
      double gain_; // gain factor (C/MeV)
  };

}
#endif

