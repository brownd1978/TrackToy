//
//  Macro to plot Mu2e muon beam info.  These are muons artifically stopped at the entrance to the Detector Solenoid
//
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
void MuStops() {
  TFile *mubeam = TFile::Open("/Users/brownd/data/sim.mu2e.MuBeamCat.MDC2020k.001201_00000000.art");
  TTree* bevents = (TTree*)mubeam->Get("Events");
  TH1F* mumom = new TH1F("mumom","Muon momentum;Momentum (MeV/c)",120,0,120.0);
  TCut bkiller("mu2e::SimParticlemv_BeamFilter__POT.obj.v_.second._stoppingCode._id==54");
  TCut bmuon("mu2e::SimParticlemv_BeamFilter__POT.obj.v_.second._pdgId==13");
  bevents->Project("mumom","mu2e::SimParticlemv_BeamFilter__POT.obj.v_.second._endMomentum.R()",bmuon&&bkiller);

  TFile *mustops = TFile::Open("/Users/brownd/data/sim.mu2e.MuminusStopsCat.MDC2020k.001201_00000000.art");
  TTree* sevents = (TTree*)mustops->Get("Events");
  TCut smuon("mu2e::SimParticlemv_TargetStopFilter__MuBeamResampler.obj.v_.second._pdgId==13");
  TCut skiller("mu2e::SimParticlemv_TargetStopFilter__MuBeamResampler.obj.v_.second._stoppingCode._id==54");
  TCut stop("mu2e::SimParticlemv_TargetStopFilter__MuBeamResampler.obj.v_.second._stoppingCode._id==32");
  TH1F* mumoms = new TH1F("mumoms","Muon momentum;Momentum (MeV/c)",120,0,120.0);
  sevents->Project("mumoms","mu2e::SimParticlemv_TargetStopFilter__MuBeamResampler.obj.v_.second._endMomentum.R()",smuon&&skiller);

  mumoms->SetLineColor(kRed);
  mumoms->SetStats(0);
  TH1F* mumomn = new TH1F("mumomn","Muon momentum;Momentum (MeV/c)",120,0,120.0);
  mumomn->Add(mumom,mumoms,5.0,-1.0);
  mumomn->SetLineColor(kBlue);
  mumomn->SetStats(0);

  TH1F* muszpos = new TH1F("muszpos","Muon stop Z position;Z (mm)",100,-420.0,420.0);
  TH2F* musxypos = new TH2F("musxypos","Muon stop XY position;X (mm);Y (mm)",40,-80.0,80,40,-80.0,80.0);
  sevents->Project("muszpos","mu2e::SimParticlemv_TargetStopFilter__MuBeamResampler.obj.v_.second._endPosition.Z()-5871",smuon&&stop);
  sevents->Project("musxypos","mu2e::SimParticlemv_TargetStopFilter__MuBeamResampler.obj.v_.second._endPosition.Y():mu2e::SimParticlemv_TargetStopFilter__MuBeamResampler.obj.v_.second._endPosition.X()+3900",smuon&&stop);
  TH1F* mustime = new TH1F("mustime","Muon stop Time;Time (ns)",100,0,500.0);
  sevents->Project("mustime","mu2e::SimParticlemv_TargetStopFilter__MuBeamResampler.obj.v_.second._endGlobalTime",smuon&&skiller);

  TCanvas* muscan = new TCanvas("MuonStops");
  muscan->Divide(2,2);
  muscan->cd(1);
  mumomn->Draw();
  mumoms->Draw("same");
  TLegend* tleg = new TLegend(0.5,0.7,0.9,0.9);
  tleg->AddEntry(mumoms,"Stopping muons","L");
  tleg->AddEntry(mumomn,"Non-stopping muons","L");
  tleg->Draw();
  muscan->cd(2);
  muszpos->Draw();
  muscan->cd(3);
  musxypos->Draw("colorz");
  muscan->cd(4);
  mustime->Draw();
}

