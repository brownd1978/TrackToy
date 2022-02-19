//
// Compare TrackToy with (TrkAna) G4 MC
//
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TH1F.h"
#include <string>
#include <iostream>

void CompareSims(const char* ttfile="Mu2eDriftAmbigCeMinusTracks.root", const char* tafile="/Users/brownd/data/TACeEMTkm.root") {
  TFile* ttf = new TFile(ttfile);
  TFile* taf = new TFile(tafile);
  TTree* trktoy = (TTree*)ttf->Get("trks");
  TTree* trkana = (TTree*)(((TDirectory*)(taf->Get("TrkAnaNeg")))->Get("trkana"));
// compute normalization
  TCut ttgood("kkstatus<=1");
  TCut tagood("de.status<=2");
  TH1F* ttmcentmom = new TH1F("ttmcentmom","MC CE momentum at tracker entrance",120,100.0,106.0);
  TH1F* tamcentmom = new TH1F("tamcentmom","MC momentum at tracker entrance",120,100.0,106.0);
  trktoy->Project("ttmcentmom","mcentmom.R()",ttgood);
  trkana->Project("tamcentmom","demcent.mom",tagood);
  double ttscale = 1.0/ttmcentmom->GetEntries();
  double tascale = 1.0/tamcentmom->GetEntries();
  cout << "Normalization for TrackToy = " << ttscale << " TrkAna " << tascale << endl;
  ttmcentmom->Scale(ttscale);
  tamcentmom->Scale(tascale);
  ttmcentmom->SetLineColor(kRed);
  tamcentmom->SetLineColor(kBlue);
  ttmcentmom->SetStats(0);
  tamcentmom->SetStats(0);

  TH1F* ttmcocost = new TH1F("ttmcocost","MC CE cos(#theta) at Origin",100,-1.0,1.0);
  TH1F* tamcocost = new TH1F("tamcocost","MC CE cos(#theta) at Origin",100,-1.0,1.0);
  trktoy->Project("ttmcocost","cos(originmom.Theta())",ttgood);
  trkana->Project("tamcocost","demc.ocosth",tagood);
  ttmcocost->Scale(ttscale);
  tamcocost->Scale(tascale);
  ttmcocost->SetLineColor(kRed);
  tamcocost->SetLineColor(kBlue);
  ttmcocost->SetStats(0);
  tamcocost->SetStats(0);

  TH1F* ttmcoposz = new TH1F("ttmcoposz","MC CE Origin Z",45,-4800,-3800);
  TH1F* tamcoposz = new TH1F("tamcoposz","MC CE Origin Z",45,-4800,-3800);
  trktoy->Project("ttmcoposz","originpos.Z()",ttgood);
  trkana->Project("tamcoposz","demc.oposz",tagood);
  ttmcoposz->Scale(ttscale);
  tamcoposz->Scale(tascale);
  ttmcoposz->SetLineColor(kRed);
  tamcoposz->SetLineColor(kBlue);
  ttmcoposz->SetStats(0);
  tamcoposz->SetStats(0);

  TCanvas* simcomp = new TCanvas("simcomp","simcomp",1000,1000);
  simcomp->Divide(2,2);
  simcomp->cd(1);
  ttmcentmom->Draw();
  tamcentmom->Draw("same");
  simcomp->cd(2);
  ttmcocost->Draw();
  tamcocost->Draw("same");
  simcomp->cd(3);
  tamcoposz->Draw();
  ttmcoposz->Draw("same");


  TH1F* ttntrkcells = new TH1F("ttntrkcells","N Track Cells",100,-0.5,99.5);
  TH1F* tantrkcells = new TH1F("tantrkcells","N Track Cells",100,-0.5,99.5);
  trktoy->Project("ttntrkcells","ncells",ttgood);
  trkana->Project("tantrkcells","de.nmatactive",tagood);
  ttntrkcells->Scale(ttscale);
  tantrkcells->Scale(tascale);
  ttntrkcells->SetLineColor(kRed);
  tantrkcells->SetLineColor(kBlue);
  ttntrkcells->SetStats(0);
  tantrkcells->SetStats(0);

  TH1F* ttntrkhits = new TH1F("ttntrkhits","N Track Active Hits",100,-0.5,99.5);
  TH1F* tantrkhits = new TH1F("tantrkhits","N Track Active Hits",100,-0.5,99.5);
  trktoy->Project("ttntrkhits","kknactive",ttgood);
  trkana->Project("tantrkhits","de.nactive",tagood);
  ttntrkhits->Scale(ttscale);
  tantrkhits->Scale(tascale);
  ttntrkhits->SetLineColor(kRed);
  tantrkhits->SetLineColor(kBlue);
  ttntrkhits->SetStats(0);
  tantrkhits->SetStats(0);

  TCanvas* trkcomp = new TCanvas("trkcomp","trkcomp",1000,1000);
  trkcomp->Divide(2,2);
  trkcomp->cd(1);
  ttntrkcells->Draw();
  tantrkcells->Draw("same");
  trkcomp->cd(2);
  ttntrkhits->Draw();
  tantrkhits->Draw("same");


}
