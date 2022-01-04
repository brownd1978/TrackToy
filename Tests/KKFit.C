#include "TFile.h"
#include "TCut.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1D.h"
#include <string>
#include <iostream>
void KKFit(const char* cefile="CeMinusTracks.root") {
  TFile* cef = new TFile(cefile);
  TTree* cet = (TTree*)cef->Get("trks");
  TH1D* momres = new TH1D("momres","KinKal Momentum Resolution at Tracker Entrance;kkmom-mcmom (MeV/c)",100,-3.0,3.0);
  TH1D* mompull = new TH1D("mompull","KinKal Momentum Pull at Tracker Entrance;(kkmom-mcmom)/kkmomerr",100,-10.0,10.0);
  TH1D* kkstatus = new TH1D("kkstatus","KinKal Fit Status",6,-0.5,5.5);
  TH1D* kkprob = new TH1D("kkprob","KinKal Fit Probability",100,0.0,1.0);
  TCut trksel ("kkstatus==0");
  cet->Project("momres","kkentmom.R()-mcentmom.R()",trksel);
  cet->Project("mompull","(kkentmom.R()-mcentmom.R())/kkentmomerr",trksel);
  cet->Project("kkstatus","kkstatus");
  cet->Project("kkprob","kkprob",trksel);
  TCanvas* kkcan = new TCanvas("kkcan","kkcan",800,800);
  kkcan->Divide(2,2);
  kkcan->cd(1);
  kkstatus->Draw();
  kkcan->cd(2);
  kkprob->Draw();
  kkcan->cd(3);
  momres->Fit("gaus");
  kkcan->cd(4);
  mompull->Fit("gaus");
}

