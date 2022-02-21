//
// utility functions for particle trajectores
//
#ifndef TrackToy_Detector_TrajUtilities_hh
#define TrackToy_Detector_TrajUtilities_hh
#include "KinKal/General/ParticleState.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include <stdexcept>
#include <stdio.h>

namespace TrackToy {
  //  Update the state of a trajectory for a change in the energy.  If this energy is still physical, this will append a new
  //  trajectory on a piecetraj at the point of the energy loss assignment and return true.  If not, it will terminate the trajectory and return false
  template<class KTRAJ> bool updateEnergy(KinKal::ParticleTrajectory<KTRAJ>& pktraj, double time, double newe, KinKal::VEC3 dmom) {
    auto const& ktraj = pktraj.nearestPiece(time);
    if(newe > pktraj.mass()) {
      // sample the momentum and position at this time
      auto momdir = (ktraj.momentum3(time) + dmom).Unit();
      auto endpos = ktraj.position3(time);
      double mass = ktraj.mass();
      // correct the momentum for the energy change
      auto newmom = sqrt(newe*newe - mass*mass)*momdir;
      // convert to particle state
      KinKal::ParticleState pstate(endpos,newmom,time,mass,ktraj.charge());
      // convert to trajectory, using the piece reference field
      KinKal::TimeRange range(time,pktraj.range().end());
      KTRAJ newtraj(pstate,ktraj.bnom(),range);
      // append this back, allowing removal
      //      std::cout << "2 appending " << range << " to range " << pktraj.range() << std::endl;
      pktraj.append(newtraj,true);
      return true;
    } else {
      // terminate the particle here
      KinKal::TimeRange range(pktraj.range().begin(),time);
      pktraj.setRange(range, true);
      return false;
    }
  }

  // extend a trajectory forwards in time till it is in tolerance with the given BField to the given z value
  template<class KTRAJ> bool extendZ(KinKal::ParticleTrajectory<KTRAJ>& pktraj, KinKal::BFieldMap const& bfield, double zmax,double tol) {
    auto tend = bfield.rangeInTolerance(pktraj.back(),pktraj.back().range().begin(),tol);
    auto pos = pktraj.position3(tend);
    while(pos.Z() < zmax && bfield.inRange(pos) && tend < pktraj.range().end() ){
      // use position and momentum at the in-tolerance end to set the next piece
      auto pstate = pktraj.state(tend);
      auto bend = bfield.fieldVect(pos);
      KinKal::TimeRange range(tend,pktraj.range().end());
      KTRAJ endtraj(pstate,bend,range);
      //        std::cout << "appending " << range << " to range " << pktraj.range() << std::endl;
      pktraj.append(endtraj);
      //        cout << "appended helix at point " << pos << " time " << range.begin() << endl;
      // update
      tend = bfield.rangeInTolerance(pktraj.back(),pktraj.back().range().begin(),tol);
      pos = pktraj.position3(tend);
    }
    return pos.Z() >= zmax;
  }

  template <class KTRAJ> void extendTraj(KinKal::BFieldMap const& bfield, KinKal::ParticleTrajectory<KTRAJ>& pktraj,double extime,double tol) {
    double tstart = pktraj.back().range().begin();
    KinKal::TimeRange range(tstart,pktraj.range().end());
    while(range.begin() < extime){
      range = KinKal::TimeRange(bfield.rangeInTolerance(pktraj.back(),range.begin(),tol),range.end());
      if(range.begin() < range.end()){
        auto pstate = pktraj.back().state(range.begin());
        auto pos = pstate.position3();
        auto bend = bfield.fieldVect(pos);
        KTRAJ endtraj(pstate,bend,range);
        //        std::cout << "appending " << range << " to range " << pktraj.range() << std::endl;
        pktraj.append(endtraj);
      } else {
        break;
      }
    }
  }

  template <class KTRAJ> double ztime(KinKal::ParticleTrajectory<KTRAJ>const& pktraj, double tstart, double zpos) {
    size_t ntries(0);
    double dz(1.0e6);
    double retval(tstart);
    do {
      ++ntries;
      auto pos = pktraj.position3(retval);
      auto vel = pktraj.velocity(retval);
      retval += (zpos-pos.Z())/vel.Z();
      pos = pktraj.position3(retval);
      dz = fabs(zpos-pos.Z());
    } while (fabs(dz)>1.0e-6 && ntries < 100);
    return retval;
  }

}
#endif
