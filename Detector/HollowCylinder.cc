#include "TrackToy/Detector/HollowCylinder.hh"
namespace TrackToy {
  std::ostream& operator <<(std::ostream& ost, HollowCylinder const& hcyl) {
    ost << "Hollow Cylinder rmin " << hcyl.rmin() << " rmax " << hcyl.rmax() << " zmin " << hcyl.zmin() << " zmax " << hcyl.zmax();
    return ost;
  }
}
