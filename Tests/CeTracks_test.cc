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
#include "KinKal/Fit/Track.hh"
#include "KinKal/Examples/ScintHit.hh" // add scint hit TODO
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
#include <fstream>
#include <vector>
#include <string>
#include <cstring>

//using namespace TrackToy;
using namespace std;
using namespace TrackToy;
using namespace KinKal;

void print_usage() {
  printf("Usage: CeTrackTest --mustopsfile s --bfieldfile s --targetrackerfile s --trackerfile s --ipafile s --endpoint f --lifetime f --tol f  --npts i --ntrks i --draw i --ttree i --minncells i --printdetail i --saveall i --cmin f --cmax f --faildetail i\n");
}

int main(int argc, char **argv) {
  using KTRAJ=LoopHelix;
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  using KKTRK = Track<KTRAJ>;
  using KKHIT = HitConstraint<KTRAJ>;
  using KKMAT = Material<KTRAJ>;
  using KKBF = BFieldEffect<KTRAJ>;
  int ntrks(-1);
  string bfile("Data/DSMapDump.dat"), mfile("MuStops.root"), targetfile("Data/Mu2eTarget.dat"), trackerfile("Data/Mu2eTracker.dat");
  string ipafile("Data/Mu2e_IPA.dat");
  string efile_my("Data/EStar_Mylar.dat"); // should come from tracker FIXME
  string sfile("Data/Schedule.txt"); // fit schedule
  double endpoint(105.0), lifetime(864.0); // these should be specified by target material FIXME
  double tol(0.01);
  double emass(0.511); //electron
  double cmin(-1.0), cmax(1.0);
  size_t npts(5000);
  bool draw(false), ttree(true), saveall(false);
  int minncells(15); // minimum # of hits
  //  double mine(90.0); // minimum energy to simulate
  // ttree variables
  TTree* cetree_(0);
  float targetde_, ipade_, trackerde_;
  float cee_, targete_, ipae_, trackere_;
  VEC3 cepos_, cemom_;
  float cet_;
  int npieces_, ntarget_, nipa_, ntrackerarcs_, ntrackercells_;
  int itrk_(0), kkstatus_, kkndof_, kknbf_, kknmat_, kknhit_, kkniter_;
  float kkchisq_, kkprob_;
  VEC3 kkentmom_, kkmidmom_, kkextmom_;
  VEC3 mcentmom_, mcmidmom_, mcextmom_;
  VEC3 kkentpos_, kkmidpos_, kkextpos_;
  VEC3 mcentpos_, mcmidpos_, mcextpos_;
  // fit parameters
  KinKal::DVEC sigmas(0.5, 0.5, 0.5, 0.5, 0.002, 0.5); // expected parameter sigmas for loop helix
  double seedsmear(10.0);
  double dwt(1.0e6);
  unsigned maxniter(10);
  Config::printLevel detail(Config::none), faildetail(Config::none);

  static struct option long_options[] = {
    {"mustopsfile",     required_argument, 0, 'm' },
    {"bfieldfile",     required_argument, 0, 'F' },
    {"targetfile",     required_argument, 0, 't' },
    {"trackerfile",     required_argument, 0, 'T' },
    {"ipafile",     required_argument, 0, 'i' },
    {"endpoint",     required_argument, 0, 'e' },
    {"tol",     required_argument, 0, 'x' },
    {"ntrks",     required_argument, 0, 'n'  },
    {"draw",     required_argument, 0, 'd'  },
    {"saveall",     required_argument, 0, 'S'  },
    {"ttree",     required_argument, 0, 'r'  },
    {"minncells",     required_argument, 0, 'M' },
    {"printdetail",     required_argument, 0, 'p' },
    {"faildetail",     required_argument, 0, 'f' },
    {"cmin",     required_argument, 0, 'c' },
    {"cmax",     required_argument, 0, 'C' },
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
      case 'N' : npts = atoi(optarg);
                 break;
      case 'n' : ntrks = atoi(optarg);
                 break;
      case 'd' : draw = atoi(optarg);
                 break;
      case 'S' : saveall = atoi(optarg);
                 break;
      case 'r' : ttree = atoi(optarg);
                 break;
      case 'M' : minncells = atoi(optarg);
                 break;
      case 'p' : detail = (Config::printLevel)atoi(optarg);
                 break;
      case 'f' : faildetail = (Config::printLevel)atoi(optarg);
                 break;
      case 'c' : cmin = atof(optarg);
                 break;
      case 'C' : cmax = atof(optarg);
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
  MatEnv::MatDBInfo matdb_(ttfinder,MatEnv::DetMaterial::moyalmean);
  BFieldMap* bfield(0);
  Config::BFCorr bfcorr(Config::variable);
  FileFinder filefinder;
  std::string fullfile = filefinder.fullFile(bfile);
  if(strncmp(bfile.c_str(),"Fixed",5) == 0){
    cout << "Building fixed BField"<< endl;
    bfield = new UniformBFieldMap(1.0);
    bfcorr = Config::nocorr;
  } else if(strncmp(bfile.c_str(),"Grad",4) == 0) {
    cout << "Building gradient BField"<< endl;
    bfield = new GradientBFieldMap(1.04,1.00, -1500,1500);
  } else {
    // setup BField
    cout << "Building BField from file " << fullfile << endl;
    bfield = new AxialBFieldMap(fullfile);
  }

  bfield->print(cout);

  // setup target
  Target target(targetfile);
  target.print(cout);
  // setup ipa
  IPA ipa(matdb_,ipafile);
  ipa.print(cout);
  // setup tracker
  Tracker tracker(matdb_,trackerfile);
  tracker.print(cout);
  EStar trackerEStar(efile_my);
  // nominal fixed field for tracking: temporary kludge FIXME
  auto bent = bfield->fieldVect(VEC3(0.0,0.0,tracker.zMin()));
  BFieldMap* trkfield = new UniformBFieldMap(bent.Z());
//  auto trkfield = bfield;
  // setup fit configuration
  Config config;
  config.dwt_ = dwt;
  config.maxniter_ = maxniter;
  config.bfcorr_ = bfcorr;
  config.tol_ = tol;
  config.plevel_ = detail;
  config.pdchi2_ = 1e6;
  // read the schedule from the file
  fullfile = filefinder.fullFile(sfile);
  std::ifstream ifs (fullfile, std::ifstream::in);
  if ( (ifs.rdstate() & std::ifstream::failbit ) != 0 ){
    std::cerr << "Error opening " << fullfile << std::endl;
    return -1;
  }
  string line;
  unsigned nmiter(0);
  while (getline(ifs,line)){
    if(strncmp(line.c_str(),"#",1)!=0){
      istringstream ss(line);
      MetaIterConfig mconfig(ss);
      mconfig.miter_ = nmiter++;
      config.schedule_.push_back(mconfig);
    }
  }
  cout << config << endl;
  // randoms
  TRandom3 tr_; // random number generator
  // setup spectrum
  CeMinusSpectrumParams ceparams(endpoint);
  CeMinusSpectrum cespect(ceparams);


  // histograms
  TFile cetracksfile("CeTracks.root","RECREATE");
  TH1F* ipade = new TH1F("ipde","IPA dE",100,-3.0,0.0);
  TH1F* tarde = new TH1F("tarde","Target dE;dE (MeV)",100,-3.0,0.0);
  TH1F* trkde = new TH1F("trkde","Tracker dE;dE (MeV)",100,-3.0,0.0);
  TH1F* trknc = new TH1F("trknc","Tracker N Cells;N Cells",100,0.001,100.0);
  if(ttree){
    cetree_ = new TTree("ce","ce");
    cetree_->Branch("itrk",&itrk_,"itrk/I");
    cetree_->Branch("cee",&cee_,"cee/F");
    cetree_->Branch("cepos",&cepos_);
    cetree_->Branch("cemom",&cemom_);
    cetree_->Branch("cet",&cet_,"cet/F");
    cetree_->Branch("ntarget",&ntarget_,"ntarget/I");
    cetree_->Branch("targetde",&targetde_,"targetde/F");
    cetree_->Branch("targete",&targete_,"targete/F");
    cetree_->Branch("ipae",&ipae_,"ipae/F");
    cetree_->Branch("nipa",&nipa_,"nipa/I");
    cetree_->Branch("ipade",&ipade_,"ipade/F");
    cetree_->Branch("npieces",&npieces_,"npieces/I");
    cetree_->Branch("ntrackerarcs",&ntrackerarcs_,"ntrackerarcs/I");
    cetree_->Branch("ntrackercells",&ntrackercells_,"ntrackercells/I");
    cetree_->Branch("trackere",&trackere_,"trackere/F");
    cetree_->Branch("trackerde",&trackerde_,"trackerde/F");
    cetree_->Branch("mcentmom",&mcentmom_);
    cetree_->Branch("mcmidmom",&mcmidmom_);
    cetree_->Branch("mcextmom",&mcextmom_);
    cetree_->Branch("mcentpos",&mcentpos_);
    cetree_->Branch("mcmidpos",&mcmidpos_);
    cetree_->Branch("mcextpos",&mcextpos_);
    cetree_->Branch("kkstatus",&kkstatus_,"kkstatus/I");
    cetree_->Branch("kkndof",&kkndof_,"kkndof/I");
    cetree_->Branch("kknbf",&kknbf_,"kknbf/I");
    cetree_->Branch("kknhit",&kknhit_,"kknhit/I");
    cetree_->Branch("kknmat",&kknmat_,"kknmat/I");
    cetree_->Branch("kkchisq",&kkchisq_,"kkchisq/F");
    cetree_->Branch("kkprob",&kkprob_,"kkprob/F");
    cetree_->Branch("kkentmom",&kkentmom_);
    cetree_->Branch("kkmidmom",&kkmidmom_);
    cetree_->Branch("kkextmom",&kkextmom_);
    cetree_->Branch("kkentpos",&kkentpos_);
    cetree_->Branch("kkmidpos",&kkmidpos_);
    cetree_->Branch("kkextpos",&kkextpos_);
  }

  std::vector<TPolyLine3D*> plhel;
  // loop over stops
  int icolor(kBlue);
  if(ntrks<0)ntrks = mtree->GetEntries();
  while (itrk_ < ntrks) {
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
    targetde_ = ipade_ = trackerde_ = 1.0;
    targete_ = ipae_ = -1.0;
    npieces_ = ntarget_ = nipa_ = ntrackerarcs_ = ntrackercells_ = -1;
    kkndof_ = kknbf_ = kknmat_ = kknhit_ = kkniter_ = -1;
    kkchisq_ = kkprob_ = -1.0;
    kkentmom_ = kkmidmom_ = kkextmom_ = VEC3();
    mcentmom_ = mcmidmom_ = mcextmom_ = VEC3();
    kkentpos_ = kkmidpos_ = kkextpos_ = VEC3();
    mcentpos_ = mcmidpos_ = mcextpos_ = VEC3();
    kkstatus_ = KinKal::Status::unfit;
    kkchisq_ = kkprob_ = -1;
    kkndof_ = kknbf_ = kknhit_ = kknmat_ = 0;

    // generate a random CeEndpoint momentum; for now, use the endpoint energy, could use spectrum TODO
    cee_ = cespect.params().EEnd_;
    // generate a random muon decay time; this should come from the target FIXME
    double tdecay = tr_.Exp(lifetime);
    // generate random phi and cos(theta)
    double phi = tr_.Uniform(-M_PI,M_PI);
    //    double cost = tr_.Uniform(0.6,0.9);
    double cost = tr_.Uniform(cmin,cmax);
    //    double cost = 0.7;
    double sint = sqrt(1.0-cost*cost);
    VEC4 const& pos4 = *mustoppos;
    cepos_ = pos4.Vect();
    cet_ = tdecay+pos4.T(); // add decay time to stopping time
    double mom = sqrt(cee_*cee_ - emass*emass);
    cemom_ = VEC3(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost);
    ParticleState cestate(cepos_,cemom_,cet_,emass,-1);
    TimeRange range(cestate.time(),cestate.time()+1000.0); // much longer than physical: is truncated later
    auto bstart = bfield->fieldVect(cestate.position3());
    KTRAJ lhelix(cestate,bstart,range);
    //    cout << "Initial trajectory " << lhelix << endl;
    // initialize piecetraj
    PKTRAJ mctraj(lhelix);
    TimeRanges targetinters, ipainters, trackerinters;
    // extend through the target
    if(target.extendTrajectory(*bfield,mctraj,targetinters)){
      //      cout << "Extended to target " << targetinters.size() << endl;
      //      mctraj.print(cout,2);
      targete_ = mctraj.energy(mctraj.range().end());
      targetde_ = targete_ - cee_;
      ntarget_ = targetinters.size();
      // extend through the IPA
      if(ipa.extendTrajectory(*bfield,mctraj,ipainters)){
        //        cout << "Extended to ipa " << ipainters.size() << endl;
        //        mctraj.print(cout,2);
        nipa_ = ipainters.size();
        ipae_ = mctraj.energy(mctraj.range().end());
        ipade_ = ipae_ - targete_;
        // extend  to the tracker
    // extend through the tracker to get the ranges
        extendZ(mctraj,*bfield, tracker.zMin()-1.0, tol);// append a fixed-field trajectory; this is a temporary kludge FIXME
        double trkent = ztime(mctraj,mctraj.back().range().begin(),tracker.zMin());
        auto endstate = mctraj.back().state(trkent);
        auto bnom = trkfield->fieldVect(endstate.position3());
        KTRAJ endtraj(endstate,bnom,KinKal::TimeRange(trkent,mctraj.range().end()));
        mctraj.append(endtraj,true); // need to understand why truncation is needed FIXME
        ////
        std::vector<double> htimes;
        std::vector<std::shared_ptr<Hit<KTRAJ>>> hits;
        std::vector<std::shared_ptr<ElementXing<KTRAJ>>> xings;
        double speed = mctraj.speed(mctraj.range().end());
        tracker.simulateHits(*trkfield,mctraj,hits,xings,trackerinters,htimes);
        ntrackercells_ = htimes.size();
        ntrackerarcs_ = trackerinters.size();
        npieces_ = mctraj.pieces().size();
        if(ntrackercells_ > minncells){
//          cout << "Extended to tracker " << trackerinters.size() << endl;
//          mctraj.print(cout,2);
// calcluate the estart energy loss
          double trackerpath(0.0);
          double ke = cestate.energy() - cestate.mass();
          for(auto const& inter : trackerinters) { trackerpath += speed*inter.range(); }
          trackerde_ = -100*trackerEStar.dEIonization(ke)*tracker.density()*trackerpath; // unit conversion
          trackere_ = mctraj.energy(mctraj.range().end());
          tarde->Fill(targetde_);
          ipade->Fill(ipade_);
          trkde->Fill(trackerde_);
          trknc->Fill(ntrackercells_);
          // add calo hits TODO
          // truncate the true trajectory
          mctraj.setRange(TimeRange(mctraj.range().begin(),mctraj.back().range().begin()+0.01));
          // get the true times at entrance and exit
          double mctent = ztime(mctraj,trackerinters.front().begin(),tracker.zMin());
          double mctext = ztime(mctraj,trackerinters.back().end(),tracker.zMax());
          double mctmid = ztime(mctraj,0.5*(mctent+mctext),tracker.zMid());
          // record true momentum at tracker entranc, mid, and exit
          auto entstate = mctraj.stateEstimate(mctent);
          auto midstate = mctraj.stateEstimate(mctmid);
          auto extstate = mctraj.stateEstimate(mctext);
          mcentmom_ = entstate.momentum3();
          mcmidmom_ = midstate.momentum3();
          mcextmom_ = extstate.momentum3();
          mcentpos_ = entstate.position3();
          mcmidpos_ = midstate.position3();
          mcextpos_ = extstate.position3();

          // reconstruct track from hits and material info
          auto seedtraj = mctraj.nearestPiece(0.5*(htimes.front()+htimes.back()));
          seedtraj.setRange(mctraj.range());
          for(size_t ipar=0; ipar < NParams(); ipar++){
            double perr = sigmas[ipar]*seedsmear;
            seedtraj.params().covariance()[ipar][ipar] = perr*perr;
            seedtraj.params().parameters()[ipar] += tr_.Gaus(0.0,perr);
          }
          if(config.plevel_ > Config::none) {
            //            cout << "Front traj " << mctraj.nearestPiece(htimes.front());
            //            cout << "Back traj " << mctraj.nearestPiece(htimes.back());
            //            cout << "Seed traj " << seedtraj;
          }
          KKTRK kktrk(config,*trkfield,seedtraj,hits,xings);
          // fill fit information
          auto const& fstat = kktrk.fitStatus();
          kkstatus_ = fstat.status_;
          kkchisq_ = fstat.chisq_.chisq();
          kkprob_ = fstat.chisq_.probability();
          kkndof_ = fstat.chisq_.nDOF();
          for(auto const& eff: kktrk.effects()) {
            const KKHIT* kkhit = dynamic_cast<const KKHIT*>(eff.get());
            const KKBF* kkbf = dynamic_cast<const KKBF*>(eff.get());
            const KKMAT* kkmat = dynamic_cast<const KKMAT*>(eff.get());
            if(kkhit != 0){
              kknhit_++;
            } else if(kkmat != 0){
              kknmat_++;
            } else if(kkbf != 0){
              kknbf_++;
            }
          }
          if(fstat.usable()){
            auto const& kktraj = kktrk.fitTraj();
            double kktent = ztime(kktraj,kktraj.range().begin(),tracker.zMin());
            double kktext = ztime(kktraj,kktraj.range().end(),tracker.zMax());
            double kktmid = ztime(kktraj,0.5*(kktent+kktext),tracker.zMid());
            auto entstate = kktraj.stateEstimate(kktent);
            auto midstate = kktraj.stateEstimate(kktmid);
            auto extstate = kktraj.stateEstimate(kktext);
            kkentmom_ = entstate.momentum3();
            kkmidmom_ = midstate.momentum3();
            kkextmom_ = extstate.momentum3();
            kkentpos_ = entstate.position3();
            kkmidpos_ = midstate.position3();
            kkextpos_ = extstate.position3();
          } else if(faildetail > Config::none) {
            cout << "Bad Fit event " << itrk_ << " status " << kktrk.fitStatus() << endl;
            cout << "True Traj " << mctraj << endl;
            cout << "Seed Traj " << seedtraj << endl;
            kktrk.print(cout,faildetail);
          }

          if(!saveall)cetree_->Fill();
          //
          if(draw){
            plhel.push_back(new TPolyLine3D(npts));
            plhel.back()->SetLineColor(icolor++%10);
            double tstart = mctraj.range().begin();
            double ts = mctraj.range().range()/(npts-1);
            VEC3 ppos;
            for(unsigned ipt=0;ipt<npts;ipt++){
              double t = tstart + ipt*ts;
              ppos = mctraj.position3(t);
              plhel.back()->SetPoint(ipt,ppos.X(),ppos.Y(),ppos.Z());
            }
          }
        } // particle hits the tracker
      } // particle stops in IPA
    } // particle exits the target going downstream
    if(saveall)cetree_->Fill();
  }
  // Canvas of basic parameters
  TCanvas* ctrkcan = new TCanvas("CeTrack");
  ctrkcan->Divide(2,2);
  ctrkcan->cd(1);
  tarde->Draw();
  ctrkcan->cd(2);
  ipade->Draw();
  ctrkcan->cd(3);
  trkde->Draw();
  ctrkcan->cd(4);
  trknc->Draw();

  ctrkcan->Write();

  if(draw){
    TCanvas* cetcan = new TCanvas("CeTracks");
    TTUBE* ttarget= new TTUBE("ttarget","ttarget","void",target.cylinder().rmin(),target.cylinder().rmax(),target.cylinder().zhalf());
    ttarget->SetLineColor(kBlack);
    ttarget->SetLineWidth(4);
    ttarget->SetFillColorAlpha(kBlack, 0.5);
    ttarget->Draw();
    TTUBE* ttracker= new TTUBE("ttracker","ttracker","void",tracker.cylinder().rmin(),tracker.cylinder().rmax(),tracker.cylinder().zhalf());
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
    rulers->SetAxisRange(target.cylinder().zmin(),tracker.cylinder().zmax(),"Z");
    rulers->Draw();
    cetcan->Write();
  }

  cetracksfile.Write();
  cetracksfile.Close();

  // now draw
  mustopsfile->Close();
  return 0;
}

