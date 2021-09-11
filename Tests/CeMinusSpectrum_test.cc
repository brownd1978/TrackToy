//
//  Test program for CeMinus Spectrum
//
#include "TrackToy/Spectra/CeMinusSpectrum.hh"
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

using namespace TrackToy;
using namespace std;

void print_usage() {
  printf("Usage: CeMinusSpectrum --endpoint f --emin f --emax f --nsample i\n");
}

int main(int argc, char **argv) {
  double emin(0.0), emax(105.0);
  double endpoint(105.0); // this should be specified by material FIXME
  size_t nsample(100000);

  static struct option long_options[] = {
    {"endpoint",     required_argument, 0, 'e' },
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
      case 'e' : endpoint = atof(optarg);
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
  TFile cefile("CeMinusSpectrum.root","RECREATE");
  // create the spectrum object
  CeMinusSpectrumParams ceparams(endpoint);
  CeMinusSpectrum cespect(ceparams);
  double inte = cespect.integral(emin,emax,10000);
  cout << "CeMinus Spectrum integral over range " << emin << " to " << emax << " = " << inte << endl;
// sample and draw the spectrum
  std::vector<double> xpts, ypts;
  xpts.reserve(nsample);
  ypts.reserve(nsample);
  double estep = (emax-emin)/(nsample-1);
  double rmax(0.0);
  for(size_t isample=0;isample<nsample;++isample){
    double energy = emin+ isample*estep;
    xpts.push_back(energy);
    double rate = cespect.rate(energy);
    ypts.push_back(rate);
    rmax = std::max(rate,rmax);
//    cout << " energy " << energy << " rate " << rate << endl;
  }
  TGraph *spectg = new TGraph(nsample,xpts.data(),ypts.data());
  spectg->SetTitle("CeMinusSpectrum;Energy (MeV);Rate (1/MeV)");
  spectg->SetLineColor(kRed);

  // plot the graph
  TCanvas* cecan = new TCanvas("cecan","CeMinusSpectrum",1000,1000);
  spectg->Draw("AL");
  // save the canvas
  cecan->Write();
  cefile.Write();
  cefile.Close();

  return 0;
}

