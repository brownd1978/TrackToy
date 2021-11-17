//
//  Test program for Spectrum
//
#include "TrackToy/Spectra/DIOSpectrum.hh"
#include "TrackToy/Spectra/CeMinusSpectrum.hh"
#include "TrackToy/Spectra/CeEndpoint.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TRandom3.h"
#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <string>

using namespace TrackToy;
using namespace std;

void print_usage() {
  printf("Usage: DIOSpectrum --file s --endpoint f --emin f --emax f --nsample i --ntest i --spectrum i\n");
}

int main(int argc, char **argv) {
  double emin(-1.0), emax(1.0e10);
  size_t nsample(100000), ntest(100000);
  double endpoint(105.0); // this should be specified by material FIXME
  string file("Data/DIOAl_fine.dat");
  int stype(Spectrum::DIO);

  static struct option long_options[] = {
    {"file",     required_argument, 0, 'f' },
    {"endpoint",     required_argument, 0, 'e' },
    {"emin",     required_argument, 0, 'm' },
    {"emax",     required_argument, 0, 'M'  },
    {"nsample",     required_argument, 0, 'N'  },
    {"ntest",     required_argument, 0, 'n'  },
    {"spectrum",     required_argument, 0, 's'  },
    {NULL, 0,0,0}
  };
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'f' : file = string(optarg);
                 break;
      case 'e' : endpoint = atof(optarg);
                 break;
      case 'm' : emin = atof(optarg);
                 break;
      case 'M' : emax = atof(optarg);
                 break;
      case 'N' : nsample = atoi(optarg);
                 break;
      case 'n' : ntest = atoi(optarg);
                 break;
      case 's' : stype = atoi(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  // setup histogramTFile diofile("DIOSpectrum.root","RECREATE");
  TFile diofile("Spectrum.root","RECREATE");
  CeMinusSpectrumParams ceparams(endpoint);
  Spectrum* spect(0);
  string title;
  switch (stype) {
    case Spectrum::DIO : default:
      spect = new DIOSpectrum(file.c_str(),emin,emax);
      title = string("DIO Spectrum;Energy (MeV);Rate (1/MeV)");
      cout << "Testing DIO spectrum" << endl;
      break;
    case Spectrum::CeMinus :
      spect = new CeMinusSpectrum(ceparams);
      title = string("CeMinus Spectrum;Energy (MeV);Rate (1/MeV)");
      cout << "Testing CeMinus spectrum" << endl;
      break;
    case Spectrum::CeEndpoint :
      spect = new CeEndpoint(endpoint);
      title = string("CeEndpoint Spectrum;Energy (MeV);Rate (1/MeV)");
      cout << "Testing CeEndpoint spectrum" << endl;
      break;
  }

  double inte = spect->integral(spect->eMin(),spect->eMax(),100000);
  cout << "Spectrum integral over range " << spect->eMin() << " to " << spect->eMax() << " = " << inte << endl;
  // sample and draw the spectrum
  std::vector<double> xpts, ypts;
  xpts.reserve(nsample);
  ypts.reserve(nsample);
  double estep = (spect->eMax()-spect->eMin())/(nsample-1);
  double rmax(0.0);
  for(size_t isample=0;isample<nsample;++isample){
    double energy = spect->eMin()+ isample*estep;
    xpts.push_back(energy);
    double rate = spect->rate(energy);
    ypts.push_back(rate);
    rmax = std::max(rate,rmax);
    //    cout << " energy " << energy << " rate " << rate << endl;
  }
  TGraph *spectg = new TGraph(nsample,xpts.data(),ypts.data());
  spectg->SetTitle(title.c_str());
  spectg->SetLineColor(kRed);
  // fill sample histogram

  TH1D* specthist = new TH1D("Spectrum",title.c_str(),1000,spect->eMin(),spect->eMax());
  specthist->SetStats(0);
  TRandom3 tr_; // random number generator
  double weight(specthist->GetNbinsX()/( spect->normalization()*ntest*(spect->eMax()-spect->eMin())));
  for(unsigned itest=0; itest<ntest; ++itest){
    double prob = tr_.Uniform();
    double energy = spect->sample(prob);
    specthist->Fill(energy,weight);
  }
  // plot the graph
  TCanvas* spectcan = new TCanvas("spectcan","Spectrum",1000,1000);
  spectcan->Divide(1,1);
  spectcan->cd(1);
  specthist->Draw();
  spectg->Draw("same");
  TLegend* leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->AddEntry(spectg,"Calculation","l");
  leg->AddEntry(specthist,"Sampling","l");
  leg->Draw();
  // save the canvas
  spectcan->Write();
  diofile.Write();
  diofile.Close();

  return 0;
}

