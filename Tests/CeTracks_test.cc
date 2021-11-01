//
//  Test Ce tracks from stopped muons
//
#include "KinKal/General/AxialBFieldMap.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/MatEnv/SimpleFileFinder.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "TrackToy/General/FileFinder.hh"
#include "TrackToy/General/TrajUtilities.hh"
#include "TrackToy/Detector/HollowCylinder.hh"
#include "TrackToy/Detector/Target.hh"
#include "TrackToy/Detector/IPA.hh"
#include "TrackToy/Detector/Tracker.hh"
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
  printf("Usage: CeTrackTest --mustopsfile s --bfieldfile s --targetrackerfile s --trackerfile s --ipafile s --endpoint f --lifetime f --tol f  --tstep f --npts i --ntrks i --draw i --ttree i --minncells i \n");
}

int main(int argc, char **argv) {
  using KTRAJ=LoopHelix;
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  int ntrks(-1);
  string bfile("Data/DSMapDump.dat"), mfile("MuStops.root"), targetfile("Data/Mu2eTarget.dat"), trackerfile("Data/Mu2eTracker.dat");
  string ipafile("Data/Mu2e_IPA.dat");
  string efile_my("Data/EStar_Mylar.dat"); // should come from tracker FIXME
  double endpoint(105.0), lifetime(864.0); // these should be specified by target material FIXME
  double tstep(0.01), tol(1e-3);
  double emass(0.511); //electron
  size_t npts(5000);
  bool draw(false), ttree(true);
  int minncells(10); // minimum # of hits
  //  double mine(90.0); // minimum energy to simulate
  // ttree variables
  TTree* cetree_;
  float cee_, targetde_, ipade_, trackerde_;
  float targete_, ipae_;
  VEC3 cepos_, cemom_;
  float cet_;
  int npieces_, ntarget_, nipa_, ntrackerarcs_, ntrackercells_;
  int itrk_(0);

  static struct option long_options[] = {
    {"mustopsfile",     required_argument, 0, 'm' },
    {"bfieldfile",     required_argument, 0, 'F' },
    {"targetfile",     required_argument, 0, 't' },
    {"trackerfile",     required_argument, 0, 'T' },
    {"ipafile",     required_argument, 0, 'i' },
    {"endpoint",     required_argument, 0, 'e' },
    {"tol",     required_argument, 0, 'x' },
    {"tstep",     required_argument, 0, 's'  },
    {"ntrks",     required_argument, 0, 'n'  },
    {"draw",     required_argument, 0, 'd'  },
    {"ttree",     required_argument, 0, 'r'  },
    {"minncells",     required_argument, 0, 'M' },
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
      case 'i' : ipafile = string(optarg);
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
      case 'r' : ttree = atoi(optarg);
                 break;
      case 'M' : minncells = atoi(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  // open the input muonstops file
  TFile* mustopsfile = TFile::Open(mfile.c_str(),"READ");
  if(!mustopsfile){
    cout << "MuStop file " << mfile << " not found: did you forget to run MuStops_test?   terminating" << endl;
    return 1;
  }
  // find the TTree in the pfile
  TTreeReader reader("MuStops",mustopsfile);
  TTreeReaderValue<VEC4> mustoppos(reader, "Pos");
  TTree* mtree = (TTree*)mustopsfile->Get("MuStops");
  cout << "MuStops TTree from file " << mfile << " has " << mtree->GetEntries() << " Entries" << endl;
  // build the materials database
  MatEnv::SimpleFileFinder ttfinder(std::string("TRACKTOY_SOURCE_DIR"),std::string("/Data/"));
  cout << "Using Materials file " << ttfinder.matMtrDictionaryFileName() << endl;
  MatEnv::MatDBInfo matdb_(ttfinder,MatEnv::DetMaterial::mpv); // use most probable value definition of eloss
  // setup BField
  FileFinder filefinder;
  std::string fullfile = filefinder.fullFile(bfile);
  AxialBFieldMap axfield(fullfile);
  cout << "axial field from file " << fullfile << " is between " << axfield.zMin() << " and " << axfield.zMax() << " with " << axfield.field().size()
    << " field values from "  << axfield.field().front() << " to "  << axfield.field().back() << endl;
  // setup target
  Target target(targetfile);
  auto const& tgtcyl = target.cylinder();

  cout << "target of material " << target.material() << " density " << target.density() << " with Z between " << tgtcyl.zmin() << " and " << tgtcyl.zmax() << " rmin " << tgtcyl.rmin() << " rmax " << tgtcyl.rmax() << endl;
  // setup ipa
  IPA ipa(matdb_,ipafile);
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
  TH1F* nipa = new TH1F("nipa","N IPA Intersections",50,-0.5,49.5);
  TH1F* ipade = new TH1F("ipde","IPA Intersection mean dE",100,0.0,0.3);
  TH1F* ipader = new TH1F("ipder","IPA Intersection sample dE",100,0.0,0.3);
  TH1F* ipades = new TH1F("ipdes","IPA Intersections Sum dE",100,0.0,0.5);
  if(ttree){
    cetree_ = new TTree("ce","ce");
    cetree_->Branch("itrk",&itrk_,"itrk/I");
    cetree_->Branch("cee",&cee_,"cee/F");
    cetree_->Branch("cepos",&cepos_);
    cetree_->Branch("cemom",&cemom_);
    cetree_->Branch("cet",&cet_,"cet/F");
    cetree_->Branch("ntarget",&ntarget_,"ntarget/I");
    cetree_->Branch("targetde",&targetde_,"targetde/F");
    cetree_->Branch("npieces",&npieces_,"npieces/I");
    cetree_->Branch("nipa",&nipa_,"nipa/I");
    cetree_->Branch("ipade",&ipade_,"ipade/F");
    cetree_->Branch("ntrackerarcs",&ntrackerarcs_,"ntrackerarcs/I");
    cetree_->Branch("ntrackercells",&ntrackercells_,"ntrackercells/I");
    cetree_->Branch("trackerde",&trackerde_,"trackerde/F");
    cetree_->Branch("targete",&targete_,"targete/F");
    cetree_->Branch("ipae",&ipae_,"ipae/F");
  }

  std::vector<TPolyLine3D*> plhel;
  // loop over stops
  int icolor(kBlue);
  if(ntrks<0)ntrks = mtree->GetEntries();
  while (itrk_ < ntrks) { // need to loop over stops if ntrks > nstops TODO
    ++itrk_;
    if(!reader.Next()){
      reader.Restart();
      if(!reader.Next()){
        cout << "Unable to rewind file" << mustopsfile << endl;
        return -2;
      }
    }
    //    cout << "Track " << itrk_ << endl;
    // reset tree variables
    targetde_ = ipade_ = trackerde_ = 0.0;
    targete_ = ipae_ = 0.0;
    nipa_ = ntrackerarcs_ = ntrackercells_ = 0;
    // generate a random CeEndpoint momentum FIXME!
    cee_ = cespect.params().EEnd_;
    // generate a random muon decay time; this should come from the target FIXME
    double tdecay = tr_.Exp(lifetime);
    // generate random phi and cos(theta)
    double phi = tr_.Uniform(-M_PI,M_PI);
    double cost = tr_.Uniform(-1.0,1.0);
    double sint = sqrt(1.0-cost*cost);
    VEC4 const& pos4 = *mustoppos;
    cepos_ = pos4.Vect();
    cet_ = tdecay+pos4.T(); // add decay time to stopping time
    //    cout << "Mustop " << pos4 << endl;
    double mom = sqrt(cee_*cee_ - emass*emass);
    cemom_ = VEC3(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost);
    ParticleState cestate(cepos_,cemom_,cet_,emass,-1);
    //    cout << "Initial Ce position" << cestate.position3() << " mom " << cestate.momentum3() << " time " << cestate.time() << endl;
    TimeRange range(cestate.time(),cestate.time()+1000.0); // replace hard-coded limit with a physical estimate FIXME
    auto bstart = axfield.fieldVect(cestate.position3());
    //      cout << "bstart " << bstart << endl;
    KTRAJ lhelix(cestate,bstart,range);
    //    cout << "Initial trajectory " << lhelix << endl;
    // initialize piecetraj
    PKTRAJ pktraj(lhelix);
    // extend to the end of the target or exiting the BField (backwards)
    double ztgt = extendZ(pktraj,axfield, axfield.zMin(), tgtcyl.zmax(), tol);
    //    cout << "Z target extend " << ztgt << endl;
    if(ztgt > tgtcyl.zmin()){
      TimeRanges targetinters, ipainters, trackerinters;
      // compute target energy loss, and update the trajectory accordingly
      bool tgtstops = target.updateTrajectory(pktraj,targetinters);
      targete_ = pktraj.energy(pktraj.range().end());
      targetde_ = targete_ - cee_;
      //      cout << "ntargetinters " << targetinters.size() << endl;
      if(!tgtstops){
        ntarget_ = targetinters.size();
        // extend through the IPA
        extendZ(pktraj,axfield, axfield.zMin(), ipa.cyl().zmax(), tol);
        // find IPA intersections
        double tstart = pktraj.range().begin();
        if(targetinters.size() > 0) tstart = targetinters.back().end();
        // extend through ipa
        extendZ(pktraj,axfield, axfield.zMin(), ipa.cyl().zmax(), tol);
        // find intersections with ipa
        ipa.cyl().intersect(pktraj,ipainters,tstart,tstep);
        nipa_ = ipainters.size();
        nipa->Fill(nipa_);
        ipade_ = 0.0;
        for(auto const& ipainter : ipainters) {
          auto eloss = ipa.energyLoss(pktraj,ipainter);
          double de = eloss.mean();
          ipade->Fill(-de);
          double der = eloss.sample(tr_.Uniform(0.0,1.0)); // currently broken, FIXME
          ipader->Fill(-der);
          //        ipade_ += der;
          ipade_ += de;
        }
        ipades->Fill(-ipade_);
        bool ipastops(false);
        if(ipainters.size() > 0){
          double energy = pktraj.energy(tstart) + ipade_;
          ipastops = updateEnergy(pktraj,ipainters.back().end(),energy);
          tstart = ipainters.back().end();
        }
        ipae_ = pktraj.energy(tstart);
        if(! ipastops) {
          // extend through tracker
          extendZ(pktraj,axfield, axfield.zMin(), trackercyl.zmax(), tol);
          // find intersections with tracker
          trackercyl.intersect(pktraj,trackerinters,tstart,tstep);
          //	cout << "ntrackerinters " << trackerinters.size() << endl;
          // check that the particle reached the tracker
          double speed = pktraj.velocity(pktraj.range().begin()).R();// assume constant speed
          npieces_ = pktraj.pieces().size();
          ntrackercells_ = tracker.nCells(speed, trackerinters);
          //	cout << "ntrackercells " << ntrackercells_ << endl;
          if(ntrackercells_ > minncells){
            ntrackerarcs_ = trackerinters.size();
            double trackerpath(0.0);
            for (auto const& range : trackerinters) trackerpath += range.range()*speed;
            double ke = cestate.energy() - cestate.mass();
            trackerde_ = trackerEStar.dEIonization(ke)*tracker.density()*trackerpath/10.0;
            tarde->Fill(-targetde_);
            trkde->Fill(trackerde_);
            trknc->Fill(ntrackercells_);
            double targetpath(0.0);
            tarpath->Fill(targetpath);
            trkpath->Fill(trackerpath);
            trktime->Fill(fmod(pktraj.range().mid(),1695.0));
            // generate tracker hits and straw interactions TODO
            // fit tracker hits TODO
            // fill fit information TODO
            //
            if(draw){
              plhel.push_back(new TPolyLine3D(npts));
              plhel.back()->SetLineColor(icolor++%10);
              double tstart = pktraj.range().begin();
              double ts = pktraj.range().range()/(npts-1);
              KinKal::VEC3 ppos;
              for(unsigned ipt=0;ipt<npts;ipt++){
                double t = tstart + ipt*ts;
                ppos = pktraj.position3(t);
                plhel.back()->SetPoint(ipt,ppos.X(),ppos.Y(),ppos.Z());
              }
            }
            cetree_->Fill();
          } // particle hits the tracker
        } // particle stops in IPA
      } // particle stops in the target
      //	cout << "particle stops in target " << endl;
    } // particle exits the target going upstream
    //      cout << "particle exits bfield " << endl;
  }
  // Draw target
  TCanvas* ctrkcan = new TCanvas("CeTrack");
  ctrkcan->Divide(2,2);
  ctrkcan->cd(1);
  tarde->Draw();
  ctrkcan->cd(2);
  ipades->Draw();
  ctrkcan->cd(3);
  trkde->Draw();
  ctrkcan->cd(4);
  trknc->Draw();

  ctrkcan->Write();

  if(draw){
    TCanvas* cetcan = new TCanvas("CeTracks");
    TTUBE* ttarget= new TTUBE("ttarget","ttarget","void",tgtcyl.rmin(),tgtcyl.rmax(),tgtcyl.zhalf());
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

