//
//  Test program for particle propagation in a BField
//
#include "KinKal/General/AxialBFieldMap.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/PiecewiseTrajectory.hh"
#include "TFile.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <string>

//using namespace TrackToy;
using namespace std;

void print_usage() {
  printf("Usage: ParticleTest --pfile s --bfile s --nsample i\n");
}

int main(int argc, char **argv) {
  size_t nsample(100000);
  string pfile, bfile;

  static struct option long_options[] = {
    {"pfile",     required_argument, 0, 'f' },
    {"bfile",     required_argument, 0, 'F' },
    {"nsample",     required_argument, 0, 'N'  },
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
      case 'N' : nsample = atoi(optarg);
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
  TTreeReaderValue<KinKal::ParticleState> pstate(reader, "particle");
  TTree* ptree = (TTree*)td->Get("nt");
  cout << "Particle TTree has " << ptree->GetEntries() << " Entries" << endl;
  // setup BField
  std::string sourcedir = getenv("TRACKTOY_SOURCE_DIR");
  std::string fullfile = sourcedir+"/"+bfile;
  KinKal::AxialBFieldMap axfield(fullfile);
  while (reader.Next()) {
    cout << "pstate status " <<  pstate.GetSetupStatus() << endl;
    cout << "Read particle with position " << pstate->position3() << " time " << pstate->time() << " momentum " << pstate->momentum3() << endl;

    // create a loop helix from this
    auto bstart = axfield.fieldVect(pstate->position3());
    KinKal::LoopHelix lhelix(*pstate,bstart);
    cout << "Loop Helix " << lhelix << endl;

  }
  // setup histogram
  TFile ptestpfile("ParticleTest.root","RECREATE");

  //  diocan->Write();
  ptestpfile.Write();
  ptestpfile.Close();
  ppfile->Close();
  return 0;
}

