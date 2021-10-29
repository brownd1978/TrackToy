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
  //  trajectory on a piecetraj at the point of the energy loss assignment and return false.  If not, it will terminate the trajectory and return true.
  template<class KTRAJ> bool updateEnergy(KinKal::ParticleTrajectory<KTRAJ>& pktraj, double time, double newe) {
    auto const& ktraj = pktraj.nearestPiece(time);
    if(newe > pktraj.mass()) {
      // sample the momentum and position at this time
      auto dir = ktraj.direction(time);
      auto endpos = ktraj.position3(time);
      double mass = ktraj.mass();
      // correct the momentum for the energy change
      auto newmom = sqrt(newe*newe - mass*mass)*dir;
      // convert to particle state
      KinKal::ParticleState pstate(endpos,newmom,time,mass,ktraj.charge());
      // convert to trajectory, using the piece reference field
      KinKal::TimeRange range(time,pktraj.range().end());
      KTRAJ newtraj(pstate,ktraj.bnom(),range);
      // append this back, allowing removal
      //      std::cout << "2 appending " << range << " to range " << pktraj.range() << std::endl;
      pktraj.append(newtraj,true);
      return false;
    } else {
    // terminate the particle here
      KinKal::TimeRange range(pktraj.range().begin(),time);
      pktraj.setRange(range, true);
      return true;
    }
  }

  // extend a trajector to the given range for the given BField to the given z range
  template<class KTRAJ> double extendZ(KinKal::ParticleTrajectory<KTRAJ>& pktraj, KinKal::BFieldMap const& bfield, double zmin, double zmax,double tol) {
    double tstart = pktraj.pieces().back().range().begin();
    auto pos = pktraj.position3(tstart);
    KinKal::TimeRange range(tstart, pktraj.range().end());
    while(pos.Z() < zmax && pos.Z() > zmin){
      range.begin() = bfield.rangeInTolerance(pktraj.back(),range.begin(),tol);
      if(range.begin() < range.end()){
        // Predict new position and momentum at this end, making linear correction for BField effects
        auto pstate = pktraj.back().state(range.begin());
        pos = pstate.position3();
        auto bend = bfield.fieldVect(pos);
        KTRAJ endtraj(pstate,bend,range);
//        std::cout << "appending " << range << " to range " << pktraj.range() << std::endl;
        pktraj.append(endtraj);
        //        cout << "appended helix at point " << pos << " time " << range.begin() << endl;
      } else {
        //        cout << "ranged out " << endl;
        break;
      }
    }
    return pos.Z();
  }
}
#endif

