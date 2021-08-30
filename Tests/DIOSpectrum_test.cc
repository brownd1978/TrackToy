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

using namespace TrackToy;
using namespace std;

void print_usage() {
  printf("Usage: DIOSpectrumTest --endpoint --emin f --emax f \n");
}

int main(int argc, char **argv) {
  double emin(0.0), emax(105.0);
  double endpoint(105.0); // this should be specified by material FIXME

  static struct option long_options[] = {
    {"endpoint",     required_argument, 0, 'e' },
    {"emin",     required_argument, 0, 'm' },
    {"emax",     required_argument, 0, 'M'  },
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
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }

// setup histograms
  TFile cefile("DIOSpectrum.root","RECREATE");
  TH1D *specth = new TH1D("specth","DIO Spectrum;E (MeV/c);d#Sigma/dE (1/MeV)",100,emin,emax);
  // create the spectrum object; this is for aluminum only!
  DIOSpectrum cespect("Data/DIO_Al.dat");
  double inte = cespect.integral(emin,emax,100000);
  cout << "DIO Spectrum integral over range " << emin << " to " << emax << " = " << inte << endl;
// sample and draw the spectrum
  static const size_t nsamples=100000;
  std::vector<double> xpts, ypts;
  xpts.reserve(nsamples);
  ypts.reserve(nsamples);
  double estep = (emax-emin)/(nsamples-1);
  double rmax(0.0);
  for(size_t isample=0;isample<nsamples;++isample){
    double energy = emin+ isample*estep;
    xpts.push_back(energy);
    double rate = cespect.rate(energy);
    ypts.push_back(rate);
    rmax = std::max(rate,rmax);
//    cout << " energy " << energy << " rate " << rate << endl;
  }
  TGraph *spectg = new TGraph(nsamples,xpts.data(),ypts.data());
  spectg->SetTitle("DIOEpectrum;Energy (MeV);Rate (1/MeV)");
  spectg->SetLineColor(kRed);
  specth->SetMinimum(0.5*xpts[1]);
  specth->SetMaximum(1.25*rmax);

  // plot the graph
  TCanvas* diocan = new TCanvas("diocan","DIOSpectrum",1000,1000);
//  specth->Draw();
  spectg->Draw("AL");
  // save the canvas
  diocan->Write();
  cefile.Write();
  cefile.Close();

  return 0;
}
 
