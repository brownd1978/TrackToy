//
//  Test muons stopping in the target (location and rate)
//
#include "KinKal/General/AxialBFieldMap.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "TrackToy/General/MuonRange.hh"
#include "TrackToy/General/FileFinder.hh"
#include "TrackToy/Detector/Target.hh"
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
#include <fstream>
#include <vector>
#include <string>

//using namespace TrackToy;
using namespace std;
using namespace TrackToy;
using namespace KinKal;

void print_usage() {
  printf("Usage: ParticleTest --mubeam s --bfieldfile s --muonrangefile s --targetfile s --beameff f--tol f  --tstep f --nbeam i\n");
}

int main(int argc, char **argv) {
  using KTRAJ=LoopHelix;
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  int nbeam(-1);
  string muonbeam("MDC2020n_10pc"), bfile("Data/DSMapDump.dat"), rfile("Data/MuonRangeAl.dat"), tfile("Data/Mu2eTarget.dat");
  double beameff(0.0043407), tstep(0.01), tol(1e-3);
  double minmass(100.0); // select muons

  static struct option long_options[] = {
    {"mubeam",     required_argument, 0, 'f' },
    {"bfieldfile",     required_argument, 0, 'F' },
    {"muonrangefile",     required_argument, 0, 'm' },
    {"targetfile",     required_argument, 0, 'T' },
    {"beameff",     required_argument, 0, 'e' },
    {"tol",     required_argument, 0, 't' },
    {"tstep",     required_argument, 0, 's'  },
    {"nbeam",     required_argument, 0, 'n'  },
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
      case 'e' : beameff = atof(optarg);
                 break;
      case 't' : tol = atof(optarg);
                 break;
      case 's' : tstep = atof(optarg);
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
  // open the input muonbeam file of muon particles artificially stopped at the entrance to the DS; this comes from the running the Mu2e software StepPointMCDumper_module on the MuBeam(Cat) dataset produced by the beam production
  FileFinder filefinder;
  string mbfile = string("Data/") + muonbeam + string("_MuBeamCat.root");
  string fullfile = filefinder.fullFile(mbfile);
  TFile* mubeamfile = TFile::Open(fullfile.c_str(),"READ");
  cout << " mbfile " << fullfile << endl;
  string simefffile = string("Data/") + muonbeam + string("_MuBeamCat_SimEff.txt");
  fullfile = filefinder.fullFile(simefffile);
  std::ifstream tgt_stream(fullfile,std::ios_base::in);
  if(tgt_stream.fail()){
    std::string errmsg = std::string("File doesn't exist" )+ fullfile;
    throw std::invalid_argument(errmsg.c_str());
  }
  string line,name;
  unsigned nb, npot;
  getline(tgt_stream, line); // skip first line
  getline(tgt_stream, line); // skip first line
  istringstream iss(line);
  iss >> name >> nb >> npot >> beameff;
  cout << "Beam muon/POT = " << beameff << endl;
  // find the TTree in the mbfile
  TDirectory* td = (TDirectory*)mubeamfile->Get("StepPointMCDumper");
  TTreeReader reader("nt",td);
  ParticleState ps;
  cout << "pstate " << ps.momentum3() << endl;
  TTreeReaderValue<ParticleState> pstate(reader, "particle");
  TTree* ptree = (TTree*)td->Get("nt");
  cout << "MuonBeam particle TTree from file " << mbfile << " has " << ptree->GetEntries() << " Entries" << endl;
  // setup BField
  fullfile = filefinder.fullFile(bfile);
  AxialBFieldMap axfield(fullfile);
  cout << "axial field from file " << fullfile << " is between " << axfield.zMin() << " and " << axfield.zMax() << " with " << axfield.field().size()
    << " field values from "  << axfield.field().front() << " to "  << axfield.field().back() << endl;
  // setup target
  Target target(tfile);
  auto const& tgtcyl = target.cylinder();
  cout << "tgtcyl between " << tgtcyl.zmin() << " and " << tgtcyl.zmax() << " rmin " << tgtcyl.rmin() << " rmax " << tgtcyl.rmax() << endl;
  // muon range table
  MuonRange muonrange(rfile.c_str(),target.density());
  cout << " muon range file " << rfile << " has density " << muonrange.density() << " and ranges " << muonrange.rangeData().size() << endl;
  int ibeam(0);
  // create a TTree for the output
  string mustopname = muonbeam + string("MuStops.root");
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
  while (ibeam < nbeam) {
    if(!reader.Next()){
      reader.Restart();
      if(!reader.Next()){
        cout << "Unable to rewind file" << mubeamfile << endl;
        return -2;
      }
    }
    // select by mass
    if(pstate->mass() > minmass){
      ++ibeam;
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
      while(pos.Z() < tgtcyl.zmax()){
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
      double tstart = pstate->time();
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
        mumoms->Fill(pstate->momentum3().R());
        muszpos->Fill(stoppos.Z()-tgtcyl.zpos());
        musxypos->Fill(stoppos.X(),stoppos.Y());
        mustime->Fill(stoppos.T());
        mustops->Fill();
      } else {
        mumomn->Fill(pstate->momentum3().R());
      }
      ++nmu;
    }
  }
  // calculate the muon stopping efficiency
  double mueff = beameff*nstopped/nmu;
  cout << "found "<< nstopped << " stopped muons out of " << nmu << ", stopping ratio = " << (float)nstopped/nmu << " stops/POT " << mueff << endl;
  // save Muon efficiencies
  string effname = muonbeam + string("_MuonStopEff.txt");
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

