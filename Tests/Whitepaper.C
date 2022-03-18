//
// Mu2e-II whitepaper comparison plots
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

void Compare(const char* ttfile="Mu2eMCAmbig_Rmue1CeMinusTracks.root",const char* tafile="/Users/brownd/data/TA_CeEMixMDC2020n_10pc.root") {
  unsigned ntaCeEGen(28000); // # CeE events generated in TrkAna
  unsigned nttCeEGen(100000); // # CeE events generated in TrkAna
  double MuBeamCatEff(0.0043407); // muon beam sampling efficiency
  double MuMinusStopEff(0.000355694294940797); // muminus stops efficiency (including prescale)
  double MuMinusStopPrescale(1000.0); // muminus stop prescale
  double tastops_per_POT=MuBeamCatEff*MuMinusStopEff*MuMinusStopPrescale;
  double ttstops_per_POT(0.00145644); // from MuStops
  double capfrac(0.609); // muon capture fraction on Al
  cout << "TA stops/POT " << tastops_per_POT << " TT stops/POT " << ttstops_per_POT << endl;
  TFile* ttf = new TFile(ttfile);
  TFile* taf = new TFile(tafile);
  TTree* trktoy = (TTree*)ttf->Get("trks");
  TTree* trkana = (TTree*)(((TDirectory*)(taf->Get("TrkAnaNeg")))->Get("trkana"));
// compute normalization
  TCut ttgood("kkstatus<=1&&kkndof>=10&&origintime%1695>450");
  TCut tagood("de.status<=2");
  TH1F* ttmcpreloss = new TH1F("ttmcpreloss","Energy Loss Before Tracker;#Delta E (MeV);Tracks/POT/50 KeV",100,-5.0,0.0);
  TH1F* tamcpreloss = new TH1F("tamcpreloss","Energy Loss Before Tracker;#Delta E (MeV);Tracks/POT/50 KeV",100,-5.0,0.0);
  trktoy->Project("ttmcpreloss","mcentmom.R()-originmom.R()",ttgood);
  trkana->Project("tamcpreloss","demcent.mom-demc.omom",tagood);
  double ttscale = ttstops_per_POT*capfrac/nttCeEGen;
  double tascale = tastops_per_POT*capfrac/ntaCeEGen;
  cout << "Normalization for TrackToy = " << ttscale << " TrkAna " << tascale << endl;
  ttmcpreloss->Scale(ttscale);
  tamcpreloss->Scale(tascale);
  ttmcpreloss->SetLineColor(kRed);
  tamcpreloss->SetLineColor(kBlue);
  ttmcpreloss->SetStats(0);
  tamcpreloss->SetStats(0);

  TH1F* ttmctrkloss = new TH1F("ttmctrkloss","Energy Loss in Tracker;#Delta E (MeV);Tracks/POT/50 KeV",100,-5.0,0.0);
  TH1F* tamctrkloss = new TH1F("tamctrkloss","Energy Loss in Tracker;#Delta E (MeV);Tracks/POT/50 KeV",100,-5.0,0.0);
  trktoy->Project("ttmctrkloss","mcextmom.R()-mcentmom.R()",ttgood);
  trkana->Project("tamctrkloss","demcxit.mom-demcent.mom",tagood);
  ttmctrkloss->Scale(ttscale);
  tamctrkloss->Scale(tascale);
  ttmctrkloss->SetLineColor(kRed);
  tamctrkloss->SetLineColor(kBlue);
  ttmctrkloss->SetStats(0);
  tamctrkloss->SetStats(0);

  TH1F* ttmcocost = new TH1F("ttmcocost","Electron cos(#theta) at Origin;;Tracks/POT",100,-1.0,1.0);
  TH1F* tamcocost = new TH1F("tamcocost","Electron cos(#theta) at Origin;;Tracks/POT",100,-1.0,1.0);
  trktoy->Project("ttmcocost","cos(originmom.Theta())",ttgood);
  trkana->Project("tamcocost","demc.ocosth",tagood);
  ttmcocost->Scale(ttscale);
  tamcocost->Scale(tascale);
  ttmcocost->SetLineColor(kRed);
  tamcocost->SetLineColor(kBlue);
  ttmcocost->SetStats(0);
  tamcocost->SetStats(0);

  TH1F* ttmcoposz = new TH1F("ttmcoposz","Electron Origin Z;Z from Tracker Center (mm);Tracks/POT",45,-4800,-3800);
  TH1F* tamcoposz = new TH1F("tamcoposz","Electron Origin Z;Z from Tracker Center (mm);Tracks/POT",45,-4800,-3800);
  trktoy->Project("ttmcoposz","originpos.Z()",ttgood);
  trkana->Project("tamcoposz","demc.oposz",tagood);
  ttmcoposz->Scale(ttscale);
  tamcoposz->Scale(tascale);
  ttmcoposz->SetLineColor(kRed);
  tamcoposz->SetLineColor(kBlue);
  ttmcoposz->SetStats(0);
  tamcoposz->SetStats(0);

  TLegend* leg = new TLegend(0.1,0.6,0.7,0.9);
  leg->AddEntry(ttmcpreloss,"TrackToy Conversion Electron","L");
  leg->AddEntry(tamcpreloss,"Geant4 Conversion Electron","L");

  TCanvas* simcomp = new TCanvas("simcomp","simcomp",1000,1000);
  simcomp->Divide(2,2);
  simcomp->cd(1);
  ttmcpreloss->Draw();
  tamcpreloss->Draw("same");
  leg->Draw();
  simcomp->cd(2);
  ttmctrkloss->Draw();
  tamctrkloss->Draw("same");
  simcomp->cd(3);
  ttmcocost->Draw();
  tamcocost->Draw("same");
  simcomp->cd(4);
  ttmcoposz->Draw();
  tamcoposz->Draw("same");

  TH1F* ttntrkcells = new TH1F("ttntrkcells","N Track Cell;;Tracks/POTs",100,-0.5,99.5);
  TH1F* tantrkcells = new TH1F("tantrkcells","N Track Cells;;Tracks/POT",100,-0.5,99.5);
  trktoy->Project("ttntrkcells","ncells",ttgood);
  trkana->Project("tantrkcells","de.nmatactive",tagood);
  ttntrkcells->Scale(ttscale);
  tantrkcells->Scale(tascale);
  ttntrkcells->SetLineColor(kRed);
  tantrkcells->SetLineColor(kBlue);
  ttntrkcells->SetStats(0);
  tantrkcells->SetStats(0);

  TH1F* ttntrkhits = new TH1F("ttntrkhits","N Track Active Hits;;Tracks/POT",100,-0.5,99.5);
  TH1F* tantrkhits = new TH1F("tantrkhits","N Track Active Hits;;Tracks/POT",100,-0.5,99.5);
  trktoy->Project("ttntrkhits","kknactive",ttgood);
  trkana->Project("tantrkhits","de.nactive",tagood);
  ttntrkhits->Scale(ttscale);
  tantrkhits->Scale(tascale);
  ttntrkhits->SetLineColor(kRed);
  tantrkhits->SetLineColor(kBlue);
  ttntrkhits->SetStats(0);
  tantrkhits->SetStats(0);

  TH1F* tttrkmomres = new TH1F("tttrkmomres","Track Momentum Resolution;#Delta p (MeV/c);Tracks/POT/50 KeV/c",100,-2.5,2.5);
  TH1F* tatrkmomres = new TH1F("tatrkmomres","Track Momentum Resolution;#Delta p (MeV/c);Tracks/POT/50 KeV/c",100,-2.5,2.5);
  trktoy->Project("tttrkmomres","kkmidmom.R()-mcmidmom.R()",ttgood);
  trkana->Project("tatrkmomres","demid.mom-demcmid.mom",tagood);
  tttrkmomres->Scale(ttscale);
  tatrkmomres->Scale(tascale);
  tttrkmomres->SetLineColor(kRed);
  tatrkmomres->SetLineColor(kBlue);
  tttrkmomres->SetStats(0);
  tatrkmomres->SetStats(0);

  TH1F* tttrkmom = new TH1F("tttrkmom","Track Momentum;Reconstructed Momentum (MeV/c);Tracks/POT/50 KeV/c",100,101,106);
  TH1F* tatrkmom = new TH1F("tatrkmom","Track Momentum;Reconstructed Momentum (MeV/c);Tracks/POT/50 KeV/c",100,101,106);
  trktoy->Project("tttrkmom","kkmidmom.R()",ttgood);
  trkana->Project("tatrkmom","demid.mom",tagood);
  tttrkmom->Scale(ttscale);
  tatrkmom->Scale(tascale);
  tttrkmom->SetLineColor(kRed);
  tatrkmom->SetLineColor(kBlue);
  tttrkmom->SetStats(0);
  tatrkmom->SetStats(0);


  TCanvas* trkcomp = new TCanvas("trkcomp","trkcomp",1000,1000);
  trkcomp->Divide(2,2);
  trkcomp->cd(1);
  ttntrkcells->Draw();
  tantrkcells->Draw("same");
  trkcomp->cd(2);
  ttntrkhits->Draw();
  tantrkhits->Draw("same");
  trkcomp->cd(3);
  gPad->SetLogy();
  tatrkmomres->Draw();
  tttrkmomres->Draw("same");
  trkcomp->cd(4);
  tttrkmom->Draw();
  tatrkmom->Draw("same");

  TCanvas* wp = new TCanvas("wp","wp",1000,500);
  wp->Divide(2,2);
  wp->cd(1);
  ttmcpreloss->Draw();
  tamcpreloss->Draw("same");
  leg->Draw();
  wp->cd(2);
  ttmctrkloss->Draw();
  tamctrkloss->Draw("same");
  wp->cd(3);
  gPad->SetLogy();
  tatrkmomres->Draw();
  tttrkmomres->Draw("same");
  wp->cd(4);
  tttrkmom->Draw();
  tatrkmom->Draw("same");

}
