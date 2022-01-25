//
//  Compare momentum resolution
//
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1F.h"
#include <string>
#include <iostream>
void MomResComp(TTree* trks) {
  using namespace std;
  TCut goodfit("kkstatus<=1");
  TCut goodmom("fabs(kkentmom.R()-mcentmom.R())<0.1");
  TCut badmom("fabs(kkentmom.R()-mcentmom.R())>1.5");

  TH1F* gprob = new TH1F("gprob","log_{10}KK fit probability",100,-10,0.0);
  TH1F* gnactive = new TH1F("gnactive","Active hits",100,0.5,100.5);
  TH1F* gactive = new TH1F("gactive","Active hit fraction",100,0.0,1.0);
  TH1F* gnull = new TH1F("gnull","Null hit fraction",100,0.0,1.0);
  TH1F* gcell = new TH1F("gcell","hit/cell fraction",100,0.7,1.0);
  TH1F* gmomerr = new TH1F("gmomerr","KK fit momentum error",100,0.0,0.8);
  gprob->SetLineColor(kBlue);
  gnactive->SetLineColor(kBlue);
  gactive->SetLineColor(kBlue);
  gnull->SetLineColor(kBlue);
  gcell->SetLineColor(kBlue);
  gmomerr->SetLineColor(kBlue);

  gprob->SetStats(0);
  gnactive->SetStats(0);
  gactive->SetStats(0);
  gnull->SetStats(0);
  gcell->SetStats(0);
  gmomerr->SetStats(0);

  TH1F* bprob = new TH1F("bprob","log_{10}KK fit probability",100,-10,0.0);
  TH1F* bnactive = new TH1F("bnactive","Active hits",100,0.5,100.5);
  TH1F* bactive = new TH1F("bactive","Active hit fraction",100,0.0,1.0);
  TH1F* bnull = new TH1F("bnull","Null hit fraction",100,0.0,1.0);
  TH1F* bcell = new TH1F("bcell","hit/cell fraction",100,0.7,1.0);
  TH1F* bmomerr = new TH1F("bmomerr","KK fit momentum error",100,0.0,0.8);
  bprob->SetLineColor(kRed);
  bnactive->SetLineColor(kRed);
  bactive->SetLineColor(kRed);
  bnull->SetLineColor(kRed);
  bcell->SetLineColor(kRed);
  bmomerr->SetLineColor(kRed);
  bprob->SetStats(0);
  bnactive->SetStats(0);
  bactive->SetStats(0);
  bnull->SetStats(0);
  bcell->SetStats(0);
  bmomerr->SetStats(0);

  trks->Project("gprob","log10(kkprob)",goodfit&&goodmom);
  trks->Project("gnactive","kknactive",goodfit&&goodmom);
  trks->Project("gactive","kknactive/kknhit",goodfit&&goodmom);
  trks->Project("gnull","kknnull/ntrkhits",goodfit&&goodmom);
  trks->Project("gcell","ntrkhits/ncells",goodfit&&goodmom);
  trks->Project("gmomerr","kkmidmomerr",goodfit&&goodmom);

  trks->Project("bprob","log10(kkprob)",goodfit&&badmom);
  trks->Project("bnactive","kknactive",goodfit&&badmom);
  trks->Project("bactive","kknactive/kknhit",goodfit&&badmom);
  trks->Project("bnull","kknnull/ntrkhits",goodfit&&badmom);
  trks->Project("bcell","ntrkhits/ncells",goodfit&&badmom);
  trks->Project("bmomerr","kkmidmomerr",goodfit&&badmom);
  bprob->Scale(10.0);
  bnactive->Scale(10.0);
  bactive->Scale(10.0);
  bnull->Scale(10.0);
  bcell->Scale(10.0);
  bmomerr->Scale(10.0);

  TCanvas* mccan = new TCanvas("mccan","mccan",1000,600);
  mccan->Divide(3,2);
  mccan->cd(1);
  gprob->Draw();
  bprob->Draw("sameh");
  mccan->cd(2);
  gactive->Draw();
  bactive->Draw("sameh");
  mccan->cd(3);
  gnull->Draw();
  bnull->Draw("sameh");
  mccan->cd(4);
  gcell->Draw();
  bcell->Draw("sameh");
  mccan->cd(5);
  gmomerr->Draw();
  bmomerr->Draw("sameh");
  mccan->cd(6);
  gnactive->Draw();
  bnactive->Draw("sameh");
}
