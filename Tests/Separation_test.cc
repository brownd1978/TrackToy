//
//  Test separating DIO from CeMinus
//
#include "TrackToy/Spectra/DIOSpectrum.hh"
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
#include <string>

using namespace TrackToy;
using namespace std;

void print_usage() {
  printf("Usage: Separation --diofile s --emin f --emax f --nsample i --stoppedmu i --separation f\n");--file s --emin f --emax f --nsample i\n");
}

int main(int argc, char **argv) {
  double emin(101.0), emax(105.0);
  size_t nsample(1000);
  double nstopped(1e22);  // should come from targetsim TODO
  double separation(10.0); // separation to define minimal detectable Rmue.  Interpretation is algorithm dependent
  double decayfrac(0.39); // muon decay fraction, should come from material info TODO
  double capfrac = 1.0-decayfrac;
  string diofile("Data/DIOAl_fine.dat");

  static struct option long_options[] = {
    {"diofile",     required_argument, 0, 'd' },
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
      case 'd' : diofile = string(optarg);
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
  TFile rootfile("Separation.root","RECREATE");
  // create the DIO spectrum object; this is for aluminum only, TODO
  DIOSpectrum diospect(diofile.c_str());
  // create the ce spectrum object
  CeMinusSpectrumParams ceparams(endpoint);
  CeMinusSpectrum cespect(ceparams);
// sample and draw the spectrum
  std::vector<double> xpts, ceypts, , dioypts;
  xpts.reserve(nsample);
  ceypts.reserve(nsample);
  dioypts.reserve(nsample);
  double estep = (emax-emin)/(nsample-1);
  double rmax(0.0);
  for(size_t isample=0;isample<nsample;++isample){
    double energy = emin+ isample*estep;
    xpts.push_back(energy);
    double rate = cespect.rate(energy);
    ceypts.push_back(rate);
    rmax = std::max(rate,rmax);
    rate = diospect.rate(energy);
    dioypts.push_back(rate);
    rmax = std::max(rate,rmax);
  }
  TGraph *rawcespectg = new TGraph(nsample,expts.data(),ceypts.data());
  rawcespectg->SetTitle("DIO Spectrum;Energy (MeV);Rate (1/MeV)");
  rawcespectg->SetLineColor(kRed);

  TGraph *rawdiospectg = new TGraph(nsample,expts.data(),dioypts.data());
  rawdiospectg->SetTitle("DIO Spectrum;Energy (MeV);Rate (1/MeV)");
  rawdiospectg->SetLineColor(kBlue);

  // plot the graph
  TCanvas* diocan = new TCanvas("diocan","DIOSpectrum",1000,1000);
  rawcespectg->Draw("AL");
  rawdiospectg->Draw("same");
  // save the canvas
  diocan->Write();
  rootfile.Write();
  rootfile.Close();

  return 0;
}
