//
//  Test program for particle propagation in a BField
//
#include "KinKal/General/CylBFieldMap.hh"
#include "TrackToy/General/FileFinder.hh"
#include "TFile.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TTree.h"
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

using namespace TrackToy;
using namespace std;

void print_usage() {
  printf("Usage: CylBFieldTest --file s --rmax f --step f\n");
}

int main(int argc, char **argv) {
  double rmax(800.0);
  double step(5.0);
  string file("Data/BMap_2D_10mm.txt");

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
  if(file.size()==0){
    cout << "No input file specified: terminating" << endl;
    return 1;
  }
  // open the input file and parse
  FileFinder filefinder;
  std::string fullfile = filefinder.fullFile(file);
  KinKal::CylBFieldMap cylfield(fullfile);
  // cout << "Loaded CylBFieldMap!\n";
  cylfield.print();
  // setup histogram
  TFile btestfile("CylBFieldTest.root","RECREATE");
  int nzbin = (int)std::ceil((cylfield.zMax()-cylfield.zMin())/step);
  int nrbin = (int)std::ceil(rmax/step);
  TH2F* brmap = new TH2F("brmap","Br field value;Z (mm); R (mm); Br (Tesla)",
      nzbin,cylfield.zMin(),cylfield.zMax(),
      nrbin,0.0,rmax);
  TH2F* bphimap = new TH2F("bphimap","Bphi field value;Z (mm); R (mm); Bphi (Tesla)",
      nzbin,cylfield.zMin(),cylfield.zMax(),
      nrbin,0.0,rmax);
  TH2F* bzmap = new TH2F("bzmap","Bz field value;Z (mm); R (mm); Bz (Tesla)",
      nzbin,cylfield.zMin(),cylfield.zMax(),
      nrbin,0.0,rmax);
  // sample the field on a grid
  for(int iz =0; iz<nzbin; ++iz){
    double zval =cylfield.zMin()+(iz+0.5)*step;
    for(int ir = 0; ir<nrbin; ++ir){
      // cout << "In ir loop\n";
      double rval =(ir+0.5)*step;
      KinKal::VEC3 rpos(rval,0.0,zval);
      KinKal::VEC3 bvec = cylfield.fieldVect(rpos);
      brmap->SetBinContent(iz+1,ir+1,bvec.X()*cos(rpos.Phi()) + bvec.Y()*sin(rpos.Phi()));
      bphimap->SetBinContent(iz+1,ir+1,-bvec.X()*sin(rpos.Phi()) + bvec.Y()*cos(rpos.Phi()));
      bzmap->SetBinContent(iz+1,ir+1,bvec.Z());
    }
  }
  TCanvas* brcan = new TCanvas("brcan","Br field map",1000,1000);
  TCanvas* bphican = new TCanvas("bphican","Bphi field map",1000,1000);
  TCanvas* bzcan = new TCanvas("bzcan","Bz field map",1000,1000);
  brcan->cd(0);
  brmap->Draw("surf");
  bphican->cd(0);
  bphimap->Draw("surf");
  bzcan->cd(0);
  bzmap->Draw("surf");
  brcan->Write();
  bphican->Write();
  bzcan->Write();
  // bprojcan->Write();
  btestfile.Write();
  btestfile.Close();
  return 0;
}