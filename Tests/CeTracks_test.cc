//
//  Test Ce tracks from stopped muons
//
#include "KinKal/General/AxialBFieldMap.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "TrackToy/General/FileFinder.hh"
#include "TrackToy/Detector/HollowCylinder.hh"
#include "TrackToy/Detector/Tracker.hh"
#include "TrackToy/Detector/EStar.hh"
#include "TrackToy/Spectra/CeMinusSpectrum.hh"
#include "TFile.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TAxis3D.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TView.h"
#include "TTUBE.h"
#include "TBRIK.h"
#include "TNode.h"
#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <string>

//using namespace TrackToy;
using namespace std;
using namespace TrackToy;
using namespace KinKal;

void print_usage() {
  printf("Usage: CeTrackTest --mustopsfile s --bfieldfile s --targetrackerfile s --trackerfile s --endpoint f --lifetime f --tol f  --tstep f --npts i --ntrks i --draw i\n");
}

int main(int argc, char **argv) {
  using KTRAJ=LoopHelix;
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  int ntrks(-1);
  string bfile("Data/DSMapDump.dat"), mfile, targetfile("Data/Mu2eTarget.dat"), trackerfile("Data/Mu2eTracker.dat");
  string efile_al("Data/EStar_Al.dat"); // should come from target FIXME
  string efile_my("Data/EStar_Mylar.dat"); // should come from tracker FIXME
  double endpoint(105.0), lifetime(864.0); // these should be specified by target material FIXME
  double tstep(0.01), tol(1e-3);
  double emass(0.511); //electron
  size_t npts(5000);
  bool draw(false);

  static struct option long_options[] = {
    {"mustopsfile",     required_argument, 0, 'm' },
    {"bfieldfile",     required_argument, 0, 'F' },
    {"targetfile",     required_argument, 0, 't' },
    {"trackerfile",     required_argument, 0, 'T' },
    {"endpoint",     required_argument, 0, 'e' },
    {"tol",     required_argument, 0, 'x' },
    {"tstep",     required_argument, 0, 's'  },
    {"ntrks",     required_argument, 0, 'n'  },
    {"draw",     required_argument, 0, 'd'  },
    {NULL, 0,0,0}
  };
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'F' : bfile = string(optarg);
                 break;
      case 'm' : mfile = string(optarg);
                 break;
      case 't' : targetfile = string(optarg);
                 break;
      case 'T' : trackerfile = string(optarg);
                 break;
      case 'e' : endpoint = atof(optarg);
                 break;
      case 'x' : tol = atof(optarg);
                 break;
      case 's' : tstep = atof(optarg);
                 break;
      case 'N' : npts = atoi(optarg);
                 break;
      case 'n' : ntrks = atoi(optarg);
                 break;
      case 'd' : draw = atoi(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  // not sure why this is necessary...
  gSystem->Load("lib/libTests.dylib");
  if(mfile.size()==0){
    cout << "No input muonstops file specified: terminating" << endl;
    return 1;
  }
  // open the input muonstops file
  TFile* mustopsfile = TFile::Open(mfile.c_str(),"READ");
  // find the TTree in the pfile
  TTreeReader reader("MuStops",mustopsfile);
  TTreeReaderValue<VEC4> pos(reader, "Pos");
  TTree* mtree = (TTree*)mustopsfile->Get("MuStops");
  cout << "MuStops TTree from file " << mfile << " has " << mtree->GetEntries() << " Entries" << endl;
  // setup BField
  FileFinder filefinder;
  std::string fullfile = filefinder.fullFile(bfile);
  AxialBFieldMap axfield(fullfile);
  cout << "axial field from file " << fullfile << " is between " << axfield.zMin() << " and " << axfield.zMax() << " with " << axfield.field().size()
    << " field values from "  << axfield.field().front() << " to "  << axfield.field().back() << endl;
  // setup target
  HollowCylinder target(targetfile);
  EStar targetEStar(efile_al);
  cout << "target between " << target.zmin() << " and " << target.zmax() << " rmin " << target.rmin() << " rmax " << target.rmax() << endl;
  // setup tracker
  Tracker tracker(trackerfile);
  EStar trackerEStar(efile_my);
  auto const& trackercyl = tracker.cylinder();
  cout << "tracker between " << trackercyl.zmin() << " and " << trackercyl.zmax() << " rmin " << trackercyl.rmin() << " rmax " << trackercyl.rmax() << endl;
  // randoms
  TRandom3 tr_; // random number generator
  // setup spectrum
  CeMinusSpectrumParams ceparams(endpoint);
  CeMinusSpectrum cespect(ceparams);

  // histograms
  TFile cetracksfile("CeTracks.root","RECREATE");
  TH1F* tarpath = new TH1F("tarpath","Target path;Path (mm)",100,10,1000.0);
  TH1F* trkpath = new TH1F("trkpath","Tracker path;Path (mm)",100,10,5000.0);
  TH1F* trktime = new TH1F("trktime","Track time;Time (ns)",100,0,2000.0);
  TH1F* trkde = new TH1F("trkde","Tracker <dE>;<dE> (MeV)",100,0.001,2.0);
  TH1F* tarde = new TH1F("tarde","Target <dE>;<dE> (MeV)",100,0.001,5.0);
  TH1F* trknc = new TH1F("trknc","Tracker N Cells;N Cells",100,0.001,100.0);

  std::vector<TPolyLine3D*> plhel;
  // loop over stops
  unsigned itrk(0);
  int icolor(kBlue);
  while (reader.Next() && (ntrks < 0 || itrk < ntrks)) {
    ++itrk;
    // generate a random CeEndpoint momentum FIXME!
    double energy = cespect.params().eMax_;
    // generate a random muon decay time FIXME
    double tdecay = tr_.Exp(lifetime);
    // generate random phi and cos(theta)
    double phi = tr_.Uniform(-M_PI,M_PI);
    double cost = tr_.Uniform(-1.0,1.0);
    double sint = sqrt(1.0-cost);
    VEC4 pos4 = *pos;
//    cout << "Mustop " << pos4 << endl;
    VEC3 cemom(energy*sint*cos(phi),energy*sint*sin(phi),energy*cost);
    ParticleState cestate(pos4.Vect(),cemom,tdecay+pos4.T(),emass,-1);
//    cout << "Initial Ce position" << cestate.position3() << " mom " << cestate.momentum3() << " time " << cestate.time() << endl;
    TimeRange range(cestate.time(),cestate.time()+1000.0); // FIXME
    auto bstart = axfield.fieldVect(cestate.position3());
    //      cout << "bstart " << bstart << endl;
    KTRAJ lhelix(cestate,bstart,range);
    //    cout << "Initial trajectory " << lhelix << endl;
    // initialize piecetraj
    PKTRAJ ptraj(lhelix);
    auto pos = cestate.position3();
    // extend to the end of the tracker or to the start of the BField
    while(pos.Z() < trackercyl.zmax() && pos.Z() > axfield.zMin()){
      range.begin() = axfield.rangeInTolerance(ptraj.back(),range.begin(),tol);
      if(range.begin() < range.end()){
        // Predict new position and momentum at this end, making linear correction for BField effects
        auto endstate = ptraj.back().state(range.begin());
        pos = endstate.position3();
        auto bend = axfield.fieldVect(pos);
        KTRAJ endhelix(endstate,bend,range);
        ptraj.append(endhelix);
//        cout << "appended helix at point " << pos << " time " << range.begin() << endl;
      } else {
//        cout << "ranged out " << endl;
        break;
      }
    }
    ptraj.back().range().end() = axfield.rangeInTolerance(ptraj.back(),ptraj.back().range().begin(),tol);
//    cout << "Particle with " << ptraj.pieces().size() << " pieces propagating from " << ptraj.position3(ptraj.range().begin())
//      << " to " << ptraj.position3(ptraj.range().end()) << endl;
    // find intersections with tracker and target
    TimeRanges targetranges, trackerranges;
    double speed = ptraj.velocity(ptraj.range().begin()).R();// assume constant speed
    if(ptraj.position3(ptraj.range().end()).Z() > trackercyl.zmax()){
      target.intersect(ptraj,targetranges,tstep);
      trackercyl.intersect(ptraj,trackerranges,tstep);
    }
    double targetpath(0.0), trackerpath(0.0);
    for (auto const& range : targetranges) targetpath += range.range()*speed;
    for (auto const& range : trackerranges) trackerpath += range.range()*speed;
    double ke = sqrt(energy*(energy + emass));
    double detarget = targetEStar.dEIonization(ke)*target.density()*targetpath/10.0; // why ionization and not total??? FIXME
    double detracker = trackerEStar.dEIonization(ke)*trackercyl.density()*trackerpath/10.0;
    double ntrkcell = tracker.nCells(speed, trackerranges);
//    cout << "detarget " << detarget << " detracker " << detracker << endl;
    tarde->Fill(detarget);
    trkde->Fill(detracker);
    trknc->Fill(ntrkcell);

//    cout << "Found " << targetranges.size() << " target ranges, path " << targetpath
//    << " and " << trackerranges.size() << " tracker ranges, path " << trackerpath << endl;
    tarpath->Fill(targetpath);
    trkpath->Fill(trackerpath);
    trktime->Fill(fmod(ptraj.range().mid(),1695.0));
    //
    if(draw){
      plhel.push_back(new TPolyLine3D(npts));
      plhel.back()->SetLineColor(icolor++%10);
      double tstart = ptraj.range().begin();
      double ts = ptraj.range().range()/(npts-1);
      KinKal::VEC3 ppos;
      for(unsigned ipt=0;ipt<npts;ipt++){
        double t = tstart + ipt*ts;
        ppos = ptraj.position3(t);
        plhel.back()->SetPoint(ipt,ppos.X(),ppos.Y(),ppos.Z());
      }
    }
  }
  // Draw target
  TCanvas* ctrkcan = new TCanvas("CeTrack");
  ctrkcan->Divide(2,2);
  ctrkcan->cd(1);
  tarde->Draw();
  ctrkcan->cd(2);
  trkde->Draw();
  ctrkcan->cd(3);
  trknc->Draw();
  ctrkcan->cd(4);
  trktime->Draw();
  ctrkcan->Write();

  if(draw){
    TCanvas* cetcan = new TCanvas("CeTracks");
    TTUBE* ttarget= new TTUBE("ttarget","ttarget","void",target.rmin(),target.rmax(),target.zhalf());
    ttarget->SetLineColor(kBlack);
    ttarget->SetLineWidth(4);
    ttarget->SetFillColorAlpha(kBlack, 0.5);
    ttarget->Draw();
    TTUBE* ttracker= new TTUBE("ttracker","ttracker","void",trackercyl.rmin(),trackercyl.rmax(),trackercyl.zhalf());
    ttracker->SetLineColor(kRed);
    ttracker->SetLineWidth(4);
    ttracker->SetFillColorAlpha(kRed, 0.5);
    ttracker->Draw();
    TTUBE *DS  = new TTUBE("DS","DS","void",0.0,800,4000);
    DS->SetVisibility(0);
    DS->Draw();
    TNode *node1 = new TNode("NODE1","NODE1","DS");
    node1->cd();
    TNode* targetnode = new TNode("targetnode","targetnode","ttarget",0.0,0.0,-4300);
    TNode* trackernode = new TNode("trackernode","trackernode","ttracker",0.0,0.0,0);
    targetnode->Draw();
    trackernode->Draw();
    node1->cd();
    node1->Draw();
    // draw tracks
    for(auto const& ph : plhel) {
      ph->Draw();
    }
    // draw the origin and axes
    TAxis3D* rulers = new TAxis3D();
    rulers->GetXaxis()->SetAxisColor(kBlue);
    rulers->GetXaxis()->SetLabelColor(kBlue);
    rulers->GetYaxis()->SetAxisColor(kCyan);
    rulers->GetYaxis()->SetLabelColor(kCyan);
    rulers->GetZaxis()->SetAxisColor(kOrange);
    rulers->GetZaxis()->SetLabelColor(kOrange);
    rulers->SetAxisRange(axfield.zMin(),trackercyl.zmax(),"Z");
    rulers->Draw();
    cetcan->Write();
  }

  cetracksfile.Write();
  cetracksfile.Close();

  // now draw
  mustopsfile->Close();
  return 0;
}

