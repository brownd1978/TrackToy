//
//  Test muons stopping in the target (location and rate)
//
#include "KinKal/General/AxialBFieldMap.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "TrackToy/General/MuonRange.hh"
#include "TrackToy/General/FileFinder.hh"
#include "TrackToy/Detector/CylinderElem.hh"
#include "TFile.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <string>

//using namespace TrackToy;
using namespace std;
using namespace TrackToy;
using namespace KinKal;

void print_usage() {
  printf("Usage: ParticleTest --particlefile s --bfieldfile s --muonrangefile s --targetfile s --tol f  --tstep f --ntrks i\n");
}

int main(int argc, char **argv) {
  using KTRAJ=LoopHelix;
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  int ntrks(-1);
  string pfile, bfile("Data/DSMapDump.dat"), rfile("Data/MuonRangeAl.dat"), tfile("Data/Mu2eTarget.dat");
  double tstep(0.01), tol(1e-3);
  double minmass(100.0); // select muons

  static struct option long_options[] = {
    {"particlefile",     required_argument, 0, 'f' },
    {"bfieldfile",     required_argument, 0, 'F' },
    {"muonrangefile",     required_argument, 0, 'm' },
    {"tol",     required_argument, 0, 't' },
    {"tstep",     required_argument, 0, 's'  },
    {"ntrks",     required_argument, 0, 'n'  },
    {NULL, 0,0,0}
  };
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'f' : pfile = string(optarg);
		 break;
      case 'F' : bfile = string(optarg);
		 break;
      case 'm' : rfile = string(optarg);
		 break;
      case 't' : tol = atof(optarg);
		 break;
      case 's' : tstep = atof(optarg);
		 break;
      case 'n' : ntrks = atoi(optarg);
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }
  // not sure why this is necessary...
  gSystem->Load("lib/libTests.dylib");
  if(pfile.size()==0){
    cout << "No input pfile specified: terminating" << endl;
    return 1;
  }
  // open the input particle pfile
  TFile* ppfile = TFile::Open(pfile.c_str(),"READ");
  // find the TTree in the pfile
  TDirectory* td = (TDirectory*)ppfile->Get("StepPointMCDumper");
  TTreeReader reader("nt",td);
  TTreeReaderValue<ParticleState> pstate(reader, "particle");
  TTree* ptree = (TTree*)td->Get("nt");
  cout << "Particle TTree from file " << pfile << " has " << ptree->GetEntries() << " Entries" << endl;
  // setup BField
  FileFinder filefinder;
  std::string fullfile = filefinder.fullFile(bfile);
  AxialBFieldMap axfield(fullfile);
  cout << "axial field from file " << fullfile << " is between " << axfield.zMin() << " and " << axfield.zMax() << " with " << axfield.field().size()
    << " field values from "  << axfield.field().front() << " to "  << axfield.field().back() << endl;
  // setup target
  CylinderElem target(tfile);
  cout << "target between " << target.zmin() << " and " << target.zmax() << " rmin " << target.rmin() << " rmax " << target.rmax() << endl;
  // muon range table
  MuonRange muonrange(rfile.c_str(),target.density());
  cout << " muon range file " << rfile << " has density " << muonrange.density() << " and ranges " << muonrange.rangeData().size() << endl;
  unsigned itrk(0);
  // histograms
  TFile mustopfile("MuStops.root","RECREATE");
  TH1F* mumoms = new TH1F("mumoms","Muon momentum;Momentum (MeV/c)",120,0,120.0);
  mumoms->SetLineColor(kRed);
  mumoms->SetStats(0);
  TH1F* mumomn = new TH1F("mumomn","Muon momentum;Momentum (MeV/c)",120,0,120.0);
  mumomn->SetLineColor(kBlue);
  mumomn->SetStats(0);
  TH1F* muszpos = new TH1F("muszpos","Muon stop Z position;Z (mm)",100,-420.0,420.0);
  TH2F* musxypos = new TH2F("musxypos","Muon stop XY position;X (mm);Y (mm)",40,-80.0,80,40,-80.0,80.0);
  unsigned nstopped(0), nmu(0);
  while (reader.Next() && (ntrks < 0 || itrk < ntrks)) {
    // select by mass
    if(pstate->mass() > minmass){
      ++itrk;
      // find the range for this muon (in mm)
      double murange = muonrange.rangeMomentum(pstate->momentum3().R());
      //      cout << "Found muon with momementum " << pstate->momentum3().R() << " range " << murange << endl;
      // build the trajectory
      TimeRange range(pstate->time(),pstate->time()+(axfield.zMax()-pstate->position3().Z())/pstate->velocity().Z());
      //      cout << "range " << range << " bpos " << pstate->position3() << endl;
      auto bstart = axfield.fieldVect(pstate->position3());
      //      cout << "bstart " << bstart << endl;
      KTRAJ lhelix(*pstate,bstart,range);
      //    cout << "Initial trajectory " << lhelix << endl;
      // initialize piecetraj
      PKTRAJ ptraj(lhelix);
      auto pos = pstate->position3();
      // extend to the end of the target
      while(pos.Z() < target.zmax()){
	range.begin() = axfield.rangeInTolerance(ptraj.back(),range.begin(),tol);
	if(range.begin() < range.end()){
	  // Predict new position and momentum at this end, making linear correction for BField effects
	  auto endstate = ptraj.back().state(range.begin());
	  pos = endstate.position3();
	  auto bend = axfield.fieldVect(pos);
	  KTRAJ endhelix(endstate,bend,range);
	  ptraj.append(endhelix);
	  //      cout << "appended helix at point " << pos << " time " << range.begin() << endl;
	} else
	  break;
      }
      //     cout << "Particle with " << ptraj.pieces().size() << " pieces propagating from " << ptraj.position3(ptraj.range().begin())
      //       << " to " << ptraj.position3(ptraj.range().end()) << endl;
      TimeRanges tranges;
      double speed = ptraj.velocity(ptraj.range().begin()).R();// assume constant speed
      target.intersect(ptraj,tranges,tstep);
      //      cout << "Found " << tranges.size() << " Intersecting ranges, with boundaries:" << endl;
      double path(0.0);
      for (auto const& range : tranges) {
	//	auto bpos = ptraj.position3(range.begin());
	//	auto epos = ptraj.position3(range.end());
	//	cout << range << " enters (r,z) " << bpos.Rho() << "," << bpos.Z() << " exits (r,z) " << epos.Rho() << "," << epos.Z() << endl;
	path += range.range()*speed;
      }
      //            cout << "intersecting path = " << path << " with range " << murange << endl;
      if(path > murange){
	++nstopped;
	mumoms->Fill(pstate->momentum3().R());
	muszpos->Fill(pos.Z()-target.zpos());
	musxypos->Fill(pos.X(),pos.Y());
      } else {
	mumomn->Fill(pstate->momentum3().R());
      }
      ++nmu;
    }
  }
  cout << "found "<< nstopped << " stopped muons out of " << nmu << " ratio = " << (float)nstopped/nmu << endl;
  // now draw
  TCanvas* muscan = new TCanvas("MuonStops");
  muscan->Divide(2,2);
  muscan->cd(1);
  mumomn->Draw();
  mumoms->Draw("same");
  TLegend* tleg = new TLegend(0.5,0.7,0.9,0.9);
  tleg->AddEntry(mumoms,"Stopping muons","L");
  tleg->AddEntry(mumomn,"Non-stopping muons","L");
  tleg->Draw();
  muscan->cd(2);
  muszpos->Draw();
  muscan->cd(3);
  musxypos->Draw("colorz");

  muscan->Write();
  mustopfile.Write();
  ppfile->Close();
  return 0;
}

