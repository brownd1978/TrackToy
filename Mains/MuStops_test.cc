//
//  Test muons stopping in the target (location and rate)
//
#include "KinKal/General/AxialBFieldMap.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/MatEnv/SimpleFileFinder.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "TrackToy/General/MuonRange.hh"
#include "TrackToy/General/FileFinder.hh"
#include "TrackToy/Detector/Target.hh"
#include "TFile.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//using namespace TrackToy;
using namespace std;
using namespace TrackToy;
using namespace KinKal;

void print_usage() {
  printf("Usage: ParticleTest --mubeam s --bfieldfile s --muonrangefile s --targetfile s --tol f  --tstep f --nbeam i --suffix s\n");
}

int main(int argc, char **argv) {
  using KTRAJ=LoopHelix;
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  using Clock = std::chrono::high_resolution_clock;
  int nbeam(1000);
  string muonbeam("MDC2020n_10pc"), bfile("Data/DSMapDump.dat"), rfile("Data/MuonRangeAl.dat"), tfile("Data/Mu2eTarget.dat");
  string suffix;
  double beameff(0.0), tstep(0.01), tol(1e-3);
  double minmass(100.0); // select muons

  static struct option long_options[] = {
    {"mubeam",     required_argument, 0, 'f' },
    {"bfieldfile",     required_argument, 0, 'F' },
    {"muonrangefile",     required_argument, 0, 'm' },
    {"targetfile",     required_argument, 0, 'T' },
    {"tol",     required_argument, 0, 't' },
    {"tstep",     required_argument, 0, 's'  },
    {"nbeam",     required_argument, 0, 'n'  },
    {"suffix",     required_argument, 0, 'S'  },
    {NULL, 0,0,0}
  };
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'f' : muonbeam = string(optarg);
                 break;
      case 'F' : bfile = string(optarg);
                 break;
      case 'm' : rfile = string(optarg);
                 break;
      case 'T' : tfile = string(optarg);
                 break;
      case 't' : tol = atof(optarg);
                 break;
      case 's' : tstep = atof(optarg);
                 break;
      case 'S' : suffix = optarg;
                 break;
      case 'n' : nbeam = atoi(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  if(muonbeam.size()==0){
    cout << "No input muon beam file specified: terminating" << endl;
    return 1;
  }
  auto now = Clock::now();
  // open the input muonbeam file of muon particles artificially stopped at the entrance to the DS; this comes from the running the Mu2e software StepPointMCDumper_module on the MuBeam(Cat) dataset produced by the beam production
  FileFinder filefinder;
  string mbfile = string("Data/") + muonbeam + string("_MuBeamCat.root");
  string fullfile = filefinder.fullFile(mbfile);
  TFile* mubeamfile = TFile::Open(fullfile.c_str(),"READ");
  cout << " mbfile " << fullfile << endl;
  string simefffile = string("Data/") + muonbeam + string("_MuBeamCat_SimEff.txt");
  fullfile = filefinder.fullFile(simefffile);
  cout << "SimEff file " << fullfile << endl;
  std::ifstream tgt_stream(fullfile,std::ios_base::in);
  if(tgt_stream.fail()){
    std::string errmsg = std::string("File doesn't exist" )+ fullfile;
    throw std::invalid_argument(errmsg.c_str());
  }
  MatEnv::SimpleFileFinder ttfinder(std::string("TRACKTOY_SOURCE_DIR"),std::string("/Data/"));
  cout << "Using Materials file " << ttfinder.matMtrDictionaryFileName() << endl;
  MatEnv::MatDBInfo matdb_(ttfinder,MatEnv::DetMaterial::moyalmean);
  string line,name;
  long long nb, npot;
  getline(tgt_stream, line); // skip first line
  getline(tgt_stream, line); // skip first line
  istringstream iss(line);
  iss >> name >> nb >> npot >> beameff;
  cout << "Beam muon/POT = " << beameff << " for dataset " << name << endl;
  // find the TTree in the mbfile
  TDirectory* td = (TDirectory*)mubeamfile->Get("StepPointMCDumper");
  TTree* ptree = (TTree*)td->Get("nt");

  Double_t        state__fArray[6];
  Double_t        time_;
  Double_t        mass_;
  Int_t           charge_;
  // List of branches
  TBranch        *b_particle_state__fArray;   //!
  TBranch        *b_particle_time_;   //!
  TBranch        *b_particle_mass_;   //!
  TBranch        *b_particle_charge_;   //!
  ptree->SetBranchAddress("state_.fArray[6]", state__fArray, &b_particle_state__fArray);
  ptree->SetBranchAddress("time_", &time_, &b_particle_time_);
  ptree->SetBranchAddress("mass_", &mass_, &b_particle_mass_);
  ptree->SetBranchAddress("charge_", &charge_, &b_particle_charge_);
  ptree->GetEntry(0);
  SVEC6 psvec(state__fArray,6);
  ParticleState ps(psvec,time_, mass_, charge_);
  cout << "pstate " << ps.momentum3() << endl;
  cout << "MuonBeam particle TTree from file " << mbfile << " has " << ptree->GetEntries() << " Entries" << endl;
  // setup BField
  fullfile = filefinder.fullFile(bfile);
  AxialBFieldMap axfield(fullfile);
  cout << "axial field from file " << fullfile << " is between " << axfield.zMin() << " and " << axfield.zMax() << " with " << axfield.field().size()
    << " field values from "  << axfield.field().front() << " to "  << axfield.field().back() << endl;
  // setup target
  auto tdiff = Clock::now() - now;
  double nsecs = std::chrono::duration_cast<std::chrono::nanoseconds>(tdiff).count();
  unsigned seed = static_cast<unsigned>(rint(nsecs));
  //  cout << "Nsecs " << nsecs << " seed " << seed << endl;
  TRandom3 tr_(seed); // random number generator
  Target target(matdb_,tfile,tr_);
//  Target target(tfile,tr_);
  auto const& tgtcyl = target.cylinder();
  target.print(cout);
  // muon range table
  MuonRange muonrange(rfile.c_str(),target.density());
  cout << " muon range file " << rfile << " has density " << muonrange.density() << " and ranges " << muonrange.rangeData().size() << endl;
  int ibeam(0);
  // create a TTree for the output
  string mustopname = muonbeam + suffix + string("MuStops.root");
  TFile mustopfile(mustopname.c_str(),"RECREATE");
  KinKal::VEC4 stoppos;
  TTree* mustops = new TTree("MuStops","MuStops",1);
  mustops->Branch("Pos",&stoppos);
  // histograms
  TH1F* mumoms = new TH1F("mumoms","Muon momentum;Momentum (MeV/c)",120,0,120.0);
  mumoms->SetLineColor(kRed);
  mumoms->SetStats(0);
  TH1F* mumomn = new TH1F("mumomn","Muon momentum;Momentum (MeV/c)",120,0,120.0);
  mumomn->SetLineColor(kBlue);
  mumomn->SetStats(0);
  TH1F* muszpos = new TH1F("muszpos","Muon stop Z position;Z (mm)",100,-420.0,420.0);
  TH1F* mustime = new TH1F("mustime","Muon stop Time;Time (ns)",100,0,500.0);
  TH2F* musxypos = new TH2F("musxypos","Muon stop XY position;X (mm);Y (mm)",40,-80.0,80,40,-80.0,80.0);
  unsigned nstopped(0), nmu(0);
  if(nbeam < 0)nbeam = ptree->GetEntries();
  int ipos(0);
  while (ibeam < nbeam) {
    ++ipos;
    if(ipos >= ptree->GetEntries())ipos=0;
    ptree->GetEntry(ipos);
    SVEC6 psvec(state__fArray,6);
    ps = ParticleState(psvec,time_, mass_, charge_);
    // select by mass
    if(ps.mass() > minmass){
      ++ibeam;
      // find the range for this muon (in mm)
      double murange = muonrange.rangeMomentum(ps.momentum3().R());
      //      cout << "Found muon with momementum " << ps.momentum3().R() << " range " << murange << endl;
      // build the trajectory
      TimeRange range(ps.time(),ps.time()+(axfield.zMax()-ps.position3().Z())/ps.velocity().Z());
      //      cout << "range " << range << " bpos " << ps.position3() << endl;
      auto bstart = axfield.fieldVect(ps.position3());
      //      cout << "bstart " << bstart << endl;
      KTRAJ lhelix(ps,bstart,range);
      //    cout << "Initial trajectory " << lhelix << endl;
      // initialize piecetraj
      PKTRAJ ptraj(lhelix);
      auto pos = ps.position3();
      // extend to the end of the target
      while(pos.Z() < tgtcyl.zmax()){
        range = KinKal::TimeRange(axfield.rangeInTolerance(ptraj.back(),range.begin(),tol),range.end());
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
      double tstart = ps.time();
      tgtcyl.intersect(ptraj,tranges,tstart,tstep);
      //      cout << "Found " << tranges.size() << " Intersecting ranges, with boundaries:" << endl;
      double path(0.0);
      bool stopped(false);
      for (auto const& range : tranges) {
        //auto bpos = ptraj.position3(range.begin());
        //auto epos = ptraj.position3(range.end());
        //cout << range << " enters (r,z) " << bpos.Rho() << "," << bpos.Z() << " exits (r,z) " << epos.Rho() << "," << epos.Z() << endl;
        path += range.range()*speed;
        if(path > murange){
          stopped = true;
          // estimate the stopping position
          double tstop = range.end() - (path-murange)/speed;
          stoppos = ptraj.position4(tstop);
//          cout << "stop at time " << tstop << " in range " << range << " position r " << stoppos.Rho() << " Z " << stoppos.Z() << endl;
          break;
        }
      }
      //            cout << "intersecting path = " << path << " with range " << murange << endl;
      if(stopped){
        ++nstopped;
        mumoms->Fill(ps.momentum3().R());
        muszpos->Fill(stoppos.Z()-tgtcyl.zpos());
        musxypos->Fill(stoppos.X(),stoppos.Y());
        mustime->Fill(stoppos.T());
        mustops->Fill();
      } else {
        mumomn->Fill(ps.momentum3().R());
      }
      ++nmu;
    }
  }
  // calculate the muon stopping efficiency
  double mueff = beameff*nstopped/nmu;
  cout << "found "<< nstopped << " stopped muons out of " << nmu << ", stopping ratio = " << (float)nstopped/nmu << " stops/POT " << mueff << endl;
  // save Muon efficiencies
  string effname = muonbeam + suffix + string("MuonStopEff.txt");
  ofstream eff_stream(effname.c_str(),std::ios_base::out);
  eff_stream << mueff << endl;
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
  muscan->cd(4);
  mustime->Draw();

  muscan->Write();
  mustopfile.Write();
  mubeamfile->Close();
  return 0;
}

