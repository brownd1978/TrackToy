//
//  Test muons stopping in the target (location and rate)
//
#include "KinKal/General/AxialBFieldMap.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "TrackToy/General/MuonRange.hh"
#include "TrackToy/General/FileFinder.hh"
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
  printf("Usage: ParticleTest --particlefile s --bfieldfile s --muonrangefile s --tol f  --tstep f --targetzpos f --targetzhalfsize f--targetrmin f --targetrmax f --tar--targetdensity f --ntrks i\n");
}

int main(int argc, char **argv) {
  using KTRAJ=LoopHelix;
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  int ntrks(-1);
  string pfile, bfile("Data/DSMapDump.dat"), rfile("Data/MuonRangeAl.dat");
  double tzpos(-4304), tzhalf(400), trmin(21.3), trmax(75.0), tdensity(1.32e-2); // Mu2e parameters
  double tstep(0.01), tol(1e-3);
  double minmass(100.0); // select muons

  static struct option long_options[] = {
    {"particlefile",     required_argument, 0, 'f' },
    {"bfieldfile",     required_argument, 0, 'F' },
    {"muonrangefile",     required_argument, 0, 'm' },
    {"tol",     required_argument, 0, 't' },
    {"tstep",     required_argument, 0, 's'  },
    {"targetzpos", required_argument, 0, 'z'  },
    {"targetzhalfsize", required_argument, 0, 'Z'  },
    {"targetrmin", required_argument, 0, 'r'  },
    {"targetrmax", required_argument, 0, 'R'  },
    {"targetdensity",     required_argument, 0, 'd'  },
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
      case 'z' : tzpos = atof(optarg);
                 break;
      case 'Z' : tzhalf = atof(optarg);
                 break;
      case 'r' : trmin = atof(optarg);
                 break;
      case 'R' : trmax = atof(optarg);
                 break;
      case 'd' : tdensity = atof(optarg);
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
  cout << "Particle TTree has " << ptree->GetEntries() << " Entries" << endl;
  // setup BField
  FileFinder filefinder;
  std::string fullfile = filefinder.fullFile(bfile);
  AxialBFieldMap axfield(fullfile);
  cout << "axial field between " << axfield.zMin() << " and " << axfield.zMax() << " with " << axfield.field().size()
    << " field values from "  << axfield.field().front() << " to "  << axfield.field().back() << endl;
  // setup target
  double tzent = tzpos-tzhalf;
  double tzexit = tzpos+tzhalf;
  cout << "target between " << tzent << " and " << tzexit << " rmin " << trmin << " rmax " << trmax << endl;
// muon range table
  MuonRange muonrange(rfile.c_str(),tdensity);
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
      // extend to the end of the BField range
      while(pos.Z() < tzexit){
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

     // search for an intersection; start with the entrance
      double tent = ptraj.zTime(tzent);
      pos = ptraj.position3(tent);
//      cout << "particle enters at " << pos << endl;
      double speed = ptraj.velocity(tent).R();// assume constant speed in integral
      double pstep = speed*tstep;
      // integrate the path inside the material
      double path(0.0);
      double ttest = tent;
      bool stops(false);
      while(pos.Z() < tzexit && (!stops)){
      // find where the particle crosses a boundary
        double rho = pos.Rho();
        bool inside = rho > trmin && rho < trmax;
        bool crosses(false);
        while(pos.Z() < tzexit && (!crosses) && (!stops)){
          ttest += tstep;
          pos = ptraj.position3(ttest);
//          cout << "ttest " << ttest << " pos " << pos << endl;
          rho = pos.Rho();
          bool newinside = rho > trmin && rho < trmax;
          crosses = newinside != inside;
          if(crosses){
//            cout << "crosses ! Rho =" << rho << " Z " << pos.Z() << endl;
            inside = newinside;
          }
          // integrate the path
          if(inside)path += pstep;
          stops =(path > murange);
        }
      }
//      cout << "intersecting path = " << path << " out of " << speed*(ttest-tent) << " with range " << murange << endl;
      if(stops){
        ++nstopped;
        mumoms->Fill(pstate->momentum3().R());
        muszpos->Fill(pos.Z()-tzpos);
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

