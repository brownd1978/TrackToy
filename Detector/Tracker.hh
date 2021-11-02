//
// Cylindrical tracker
//
#ifndef TrackToy_Detector_Tracker_hh
#define TrackToy_Detector_Tracker_hh
#include "TrackToy/Detector/HollowCylinder.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Detector/StrawMaterial.hh"
#include <string>
#include <vector>
#include <iostream>
namespace TrackToy {
  class Tracker {
    public:
      enum CellOrientation{azimuthal=0,axial};
      Tracker(MatEnv::MatDBInfo const& matdbinfo,std::string const& tfile);
      auto const& cylinder() const { return cyl_; }
      unsigned nCells() const { return ncells_; }
      double density() const { return density_;} // gm/cm^3
      double cellDensity() const { return cellDensity_;} // N/mm
      void print(std::ostream& os) const;
      const KinKal::StrawMaterial* strawMaterial() const { return smat_; }
      // define hit times from a range
      // create hit and straw crossings from a trajectory
      template<class PKTRAJ> bool simulateHits(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj,
          std::vector<KinKal::TimeRange>& tinters, std::vector<double>& htimes) const;
    private:
      HollowCylinder cyl_; // geometric form of the tracker
      CellOrientation orientation_; // orientation of the cells
      unsigned ncells_; // number of cells
      double cellDensity_; // effective linear cell density mm^-1
      double density_; // total average density gm/mm^3
      KinKal::StrawMaterial* smat_; // straw material
  };

  template<class PKTRAJ> bool Tracker::simulateHits(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj,
      std::vector<KinKal::TimeRange>& tinters, std::vector<double>& htimes) const {
    double tstart = pktraj.back().range().begin();
    double speed = pktraj.speed(tstart);
    double tol = 1.0/speed; // 1mm accuracy
    double tstep = smat_->strawRadius()/speed;
    // extend through the tracker to get the ranges
    extendZ(pktraj,bfield, bfield.zMin(), cylinder().zmax(), tol);
    // find intersections with tracker
    cylinder().intersect(pktraj,tinters,tstart,tstep);
//    std::cout << "ninters " << tinters.size() << std::endl;
    for(auto const& tinter : tinters) {
      double clen(0.0);
      double time = tinter.begin();
      while(time < tinter.end()){
        auto vel = pktraj.velocity(time);
        if(orientation_ == azimuthal){
          auto pos = pktraj.position3(time);
          auto rdir = KinKal::VEC3(pos.X(),pos.Y(),0.0).Unit(); // radial direction
          double vr = vel.Dot(rdir); // radial component of velocity
          double vtot = sqrt(vel.Z()*vel.Z() + vr*vr);
          clen += vtot*tstep;
        } else {
          clen += vel.R()*tstep;
        }
        time += tstep;
      }
      unsigned ncells = (unsigned)rint(clen*cellDensity_);
//      std::cout << "clen " << clen << std::endl;
      double hstep = tinter.range()/(ncells+1);
      double htime = tinter.begin()+0.5*tstep;
      for(unsigned icell=0;icell<ncells;++icell){
        htimes.push_back(htime);
        htime += hstep;
      }
    }
    return true;
  }
}
#endif

