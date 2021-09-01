//
//  Test program for DIO Spectrum
//
#include "TrackToy/Spectra/DIOSpectrum.hh"
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
  printf("Usage: DIOSpectrum --file s --emin f --emax f --nsample i\n");
}

int main(int argc, char **argv) {
  double emin(0.0), emax(105.0);
  size_t nsample(100000);
  string file("Data/DIOAl_fine.dat");

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
  TFile diofile("DIOSpectrum.root","RECREATE");
  // create the spectrum object; this is for aluminum only!
  DIOSpectrum diospect(file.c_str());
  double inte = diospect.integral(emin,emax,100000);
  cout << "DIO Spectrum integral over range " << emin << " to " << emax << " = " << inte << endl;
// sample and draw the spectrum
  std::vector<double> xpts, ypts;
  xpts.reserve(nsample);
  ypts.reserve(nsample);
  double estep = (emax-emin)/(nsample-1);
  double rmax(0.0);
  for(size_t isample=0;isample<nsample;++isample){
    double energy = emin+ isample*estep;
    xpts.push_back(energy);
    double rate = diospect.rate(energy);
    ypts.push_back(rate);
    rmax = std::max(rate,rmax);
//    cout << " energy " << energy << " rate " << rate << endl;
  }
  TGraph *spectg = new TGraph(nsample,xpts.data(),ypts.data());
  spectg->SetTitle("DIO Spectrum;Energy (MeV);Rate (1/MeV)");
  spectg->SetLineColor(kRed);
  // plot the graph
  TCanvas* diocan = new TCanvas("diocan","DIOSpectrum",1000,1000);
  spectg->Draw("AL");
  // save the canvas
  diocan->Write();
  diofile.Write();
  diofile.Close();

  return 0;
}

