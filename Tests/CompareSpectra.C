#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1D.h"
#include <string>
#include <iostream>
void CompareSpectra(const char* diofile="FlatDIOTracks.root",double dioscale=1.0,const char* cefile="CeMinusTracks.root",double cescale=1.0, double momcut=103.7) {
  using namespace std;
  TFile* diof = new TFile(diofile);
  TTree* diot = (TTree*)diof->Get("trks");
  TFile* cef = new TFile(cefile);
  TTree* cet = (TTree*)cef->Get("trks");
  double momlow(103.0);
  double momhigh(106.0);
  unsigned nbins(30);
  TH1D* diomom = new TH1D("diomom","KinKal Momentum at Tracker Entrance;mom (MeV/c);Events/MeV/c",nbins,momlow,momhigh);
  TH1D* cemom = new TH1D("cemom","KinKal Momentum at Tracker Entrance;mom (MeV/c);Events/MeV/c",nbins,momlow,momhigh);
  diomom->SetStats(0);
  diomom->SetLineColor(kBlue);
  cemom->SetStats(0);
  cemom->SetLineColor(kRed);
  double timecut(700.0), momerrcut(0.3), hitfraccut(0.9);
  char cutstring[100];
  snprintf(cutstring,100,"(kkstatus==0&&fmod(kkmidt0,1695)>%3.1f&&kkmidmomerr<%3.3f&&ntrkhits/ncells>%3.3f)*%3.3f",timecut,momerrcut,hitfraccut,dioscale);// filter on fit quality, tracktime, etc
  cout << "DIO selection " << cutstring << endl;
  TCut diosel = TCut("weight")*TCut(cutstring); //weight to get physical rate
  snprintf(cutstring,100,"(kkstatus==0&&fmod(kkmidt0,1695)>%3.1f&&kkmidmomerr<%3.3f&&ntrkhits/ncells>%3.3f)*%3.3f",timecut,momerrcut,hitfraccut,cescale);// filter on fit quality, tracktime, etc
  cout << "Ce selection " << cutstring << endl;
  TCut cesel = TCut("weight")*TCut(cutstring); //weight to get physical rate

  diot->Project("diomom","kkentmom.R()",diosel);
  cet->Project("cemom","kkentmom.R()",cesel);
  double max = 3*cemom->GetBinContent(cemom->GetMaximumBin());
  cemom->SetMaximum(max);
  diomom->SetMaximum(max);
// cut and count
  double ndio(0.0), nce(0.0);
  for(int ibin=1;ibin<=nbins;ibin++){
    if(cemom->GetBinLowEdge(ibin) >= momcut){
      nce += cemom->GetBinContent(ibin);
      ndio += diomom->GetBinContent(ibin);
    }
  }
  cout << "NCe = " << nce << " NDIO = " << ndio << endl;
  TCanvas* ccan = new TCanvas("ccan","Compare",800,800);
//  ccan->SetLogy();
  diomom->SetMinimum(1e-3);
  diomom->Draw("h");
  cemom->Draw("hsame");
  TLegend* leg = new TLegend(0.4,0.7,0.9,0.9);
  leg->AddEntry(diomom,diofile,"l");
  leg->AddEntry(cemom,cefile,"l");
  leg->Draw();
  TPaveText* info = new TPaveText(0.4,0.5,0.9,0.7,"NDCNB");
  char tline[80];
  snprintf(tline,80,"t0>%3.1f",timecut);
  info->AddText(tline);
  snprintf(tline,80,"momerr<%3.3f",momerrcut);
  info->AddText(tline);
  snprintf(tline,80,"hitfrac>%3.3f",hitfraccut);
  info->AddText(tline);
  snprintf(tline,80,"KKEntMom>%3.1f, NCe=%3.1f, NDIO=%3.1f",momcut,nce,ndio);
  info->AddText(tline);
  info->Draw();
}
