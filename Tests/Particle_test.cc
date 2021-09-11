//
//  Test program for particle propagation in a BField
//
#include "KinKal/General/BFieldMap.hh"
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
  printf("Usage: ParticleTest --file s --zmin f --zmax f --bgrad f --nsample i\n");
}

int main(int argc, char **argv) {
  double zmin(-1500), zmax(1500.0);
  double bgrad (-0.03);
  size_t nsample(100000);
  string file;

  static struct option long_options[] = {
    {"file",     required_argument, 0, 'f' },
    {"zmin",     required_argument, 0, 'm' },
    {"zmax",     required_argument, 0, 'M'  },
    {"bgrad",     required_argument, 0, 'b'  },
    {"nsample",     required_argument, 0, 'N'  },
    {NULL, 0,0,0}
  };
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'f' : file = string(optarg);
		 break;
      case 'm' : zmin = atof(optarg);
		 break;
      case 'M' : zmax = atof(optarg);
		 break;
      case 'b' : bgrad = atof(optarg);
		 break;
      case 'N' : nsample = atoi(optarg);
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }
// not sure why this is necessary...
  gSystem->Load("lib/libTests.dylib");
// open the input particle file
  TFile* pfile = TFile::Open(file.c_str(),"READ");
  // find the TTree in the file
  TDirectory* td = (TDirectory*)pfile->Get("StepPointMCDumper");
  TTreeReader reader("nt",td);
  TTreeReaderValue<KinKal::ParticleState> pstate(reader, "particle");

  TTree* ptree = (TTree*)td->Get("nt");
  cout << "Particle TTree has " << ptree->GetEntries() << " Entries" << endl;
//  KinKal::ParticleState pstate;
//  TBranch *part = ptree->GetBranch("particle");
//  part->SetAddress(&pstate);
//  for(unsigned ipart=0;ipart<ptree->GetEntries(); ++ipart){
while (reader.Next()) {
 cout << "pstate status " <<  pstate.GetSetupStatus() << endl;
//part->GetEntry(ipart);
    cout << "Read particle with position " << pstate->position3() << " time " << pstate->time() << " momentum " << pstate->momentum3() << endl;
  }
// setup histogram
  TFile ptestfile("ParticleTest.root","RECREATE");

//  diocan->Write();
  ptestfile.Write();
  ptestfile.Close();
  pfile->Close();
  return 0;
}

