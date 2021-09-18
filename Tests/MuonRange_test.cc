//
//  Test program for Muon range
//
#include "TrackToy/General/MuonRange.hh"
#include "TFile.h"
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

using namespace TrackToy;
using namespace std;

void print_usage() {
  printf("Usage: MuonRange --file s --emin f --emax f --density f --nsample i\n");
}

int main(int argc, char **argv) {
  double emin(0.0), emax(200.0);
  size_t nsample(10000);
  string file("Data/MuonRangeAl.dat");

  static struct option long_options[] = {
    {"file",     required_argument, 0, 'f' },
    {"emin",     required_argument, 0, 'm' },
    {"emax",     required_argument, 0, 'M'  },
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
      case 'm' : emin = atof(optarg);
		 break;
      case 'M' : emax = atof(optarg);
		 break;
      case 'N' : nsample = atoi(optarg);
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }

// setup histograms
  TFile murangefile("MuonRange.root","RECREATE");
  // create the object
  MuonRange murange(file.c_str(),1.0); // unit density: return gm/cm^2
  cout << " muon range file " << file << " has density " << murange.density() << " and ranges " << murange.rangeData().size() << endl;
// sample and draw the range
  std::vector<double> xpts, epts, mpts;
  xpts.reserve(nsample);
  epts.reserve(nsample);
  mpts.reserve(nsample);
  double estep = (emax-emin)/(nsample-1);
  for(size_t isample=0;isample<nsample;++isample){
    double energy = emin+ isample*estep;
    xpts.push_back(energy);
    double erange = murange.rangeEnergy(energy);
    epts.push_back(erange);
    double mrange = murange.rangeMomentum(energy);
    mpts.push_back(mrange);
//    cout << " energy " << energy << " range " << range << endl;
  }
  TGraph *espectg = new TGraph(nsample,xpts.data(),epts.data());
  espectg->SetTitle("Muon Range;Energy (MeV);Range (gm/cm^{2})");
  espectg->SetLineColor(kRed);

  TGraph *mspectg = new TGraph(nsample,xpts.data(),mpts.data());
  mspectg->SetTitle("Muon Range;Momentum (MeV/c);Range (gm/cm^{2})");
  mspectg->SetLineColor(kBlue);

  // plot the graph
  TCanvas* mucan = new TCanvas("mucan","MuonRange",1000,1000);
  mucan->Divide(2,1);
  mucan->cd(1);
  gPad->SetLeftMargin(0.15);
  espectg->Draw("AL");
  mucan->cd(2);
  gPad->SetLeftMargin(0.15);
  mspectg->Draw("AL");
  // save the canvas
  mucan->Write();
  murangefile.Write();
  murangefile.Close();

  return 0;
}

