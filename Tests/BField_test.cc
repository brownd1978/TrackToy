//
//  Test program for particle propagation in a BField
//
#include "KinKal/General/AxialBFieldMap.hh"
#include "TFile.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TBranch.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

//using namespace TrackToy;
using namespace std;

void print_usage() {
  printf("Usage: BFieldTest --file s --rmax f --step f\n");
}

int main(int argc, char **argv) {
  double rmax(800.0);
  double step(5.0);
  string file;

  static struct option long_options[] = {
    {"file",     required_argument, 0, 'f' },
    {"rmax",     required_argument, 0, 'r' },
    {"step",     required_argument, 0, 's' },
    {NULL, 0,0,0}
  };
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'f' : file = string(optarg);
		 break;
      case 'r' : rmax = atof(optarg);
		 break;
      case 's' : step = atof(optarg);
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }
// not sure why this is necessary...
  gSystem->Load("lib/libTests.dylib");
  if(file.size()==0){
    cout << "No input file specified: terminating" << endl;
    return 1;
  }
  // open the input file and parse
  std::string sourcedir = getenv("TRACKTOY_SOURCE_DIR");
  std::string fullfile = sourcedir+"/"+file;
  // read the file
  std::ifstream spectrum_stream(fullfile);
  std::string line;
  bool first(true);
  static std::string comment("#");
  double zmax, bz, zmin, zold;
  vector<double> axial;
  axial.reserve(1000);
  while (std::getline(spectrum_stream, line)) {
    // skip comments and blank lines
    if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
      std::istringstream iss(line);
//      cout << "read line " << line << endl;
      if (iss >> zmax >> bz ) {
//	cout << "read z " << zmax << " b " << bz << endl;
	if(first){
	  first = false;
	  zmin = zmax;
	}
	if(axial.size()>0){
	  // check uniformity
	  double dz = zmax - (zmin + (zmax-zold)*axial.size());
	  if(fabs(dz) > 1e-5) cout << "Error: field spec isn't uniform!" << dz << " " << zmax << " " << axial.size() << endl;
	}
	zold = zmax;
	axial.push_back(bz);
      }
    }
  }
  cout << "Read axial field between " << zmin << " and " << zmax << " with " << axial.size() << " values: " << endl;
//  for(auto bz : axial) cout << "Bz =" << bz << endl;

  KinKal::AxialBFieldMap axfield(zmin,zmax,axial);
  // setup histogram
  TFile btestfile("BFieldTest.root","RECREATE");
  int nzbin = (int)std::ceil((zmax-zmin)/step);
  int nrbin = (int)std::ceil(rmax/step);
  TH1F* bzproj = new TH1F("bzproj","Bz field value;Z (mm); Bz (Tesla)", nzbin,zmin,zmax);
  TH1F* brproj = new TH1F("brproj","Br field value at rmax;Z (mm); Br (Tesla)", nzbin,zmin,zmax);
  TH2F* brmap = new TH2F("brmap","Br field value;Z (mm); R (mm); Br (Tesla)",
      nzbin,zmin,zmax,
      nrbin,0.0,rmax);
  // sample the field on a grid
  for(int iz =0; iz<nzbin; ++iz){
    double zval =zmin+(iz+0.5)*step;
    KinKal::VEC3 zpos(rmax,0.0,zval);
    KinKal::VEC3 bvec = axfield.fieldVect(zpos);
    bzproj->SetBinContent(iz+1,bvec.Z());
    brproj->SetBinContent(iz+1,bvec.Rho());
    for(int ir = 0; ir<nrbin; ++ir){
      double rval =(ir+0.5)*step;
      KinKal::VEC3 rpos(rval,0.0,zval);
      bvec = axfield.fieldVect(rpos);
      brmap->SetBinContent(iz+1,ir+1,bvec.Rho());
    }
  }
  TCanvas* brcan = new TCanvas("brcan","Br field map",1000,1000);
  TCanvas* bprojcan = new TCanvas("bprojcan","Bz field map",1000,1000);
  brcan->cd(0);
  brmap->Draw("surf");
  bprojcan->Divide(1,2);
  bprojcan->cd(1);
  bzproj->Draw("poly");
  bprojcan->cd(2);
  brproj->Draw("poly");
  brcan->Write();
  bprojcan->Write();
  btestfile.Write();
  btestfile.Close();
  return 0;
}

