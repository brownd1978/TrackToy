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

  template<class KTRAJ> double timeStep(KinKal::ParticleTrajectory<KTRAJ>const& pktraj, double zmin, double zmax, double tstart, double tstep) {
    double ttest = tstart;
    auto pos = pktraj.position3(ttest);
    double zpos = 0.5*(zmin+zmax);
    if(pos.Z() > zmin && pos.Z() < zmax){
      // take small steps if we're in the right z range
      ttest += tstep;
    } else {
      // step through the pieces till we are going in the right direction
      auto vel = pktraj.velocity(ttest);
      while ( (pos.Z()-zpos)*vel.z() > 0.0 && ttest < pktraj.range().end() ){
	auto const& piece = pktraj.nearestPiece(ttest);
	ttest = piece.range().end() + tstep;
	pos = pktraj.position3(ttest);
	vel = pktraj.velocity(ttest);
      }
      // step through the pieces till we're in the right z range
      while( ttest < pktraj.range().end() &&
	  ((vel.Z() > 0.0 && pos.Z() < zmin) ||
	   (vel.Z() < 0.0 && pos.Z() > zmax)) ){
	auto const& piece = pktraj.nearestPiece(ttest);
	ttest = piece.range().end() + tstep;
	pos = pktraj.position3(ttest);
      }
    }
    return ttest;
  }
}
#endif
