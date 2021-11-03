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

//using namespace TrackToy;
using namespace std;
using namespace TrackToy;
using namespace KinKal;

void print_usage() {
  printf("Usage: CeTrackTest --mustopsfile s --bfieldfile s --targetrackerfile s --trackerfile s --ipafile s --endpoint f --lifetime f --tol f  --tstep f --npts i --ntrks i --draw i --ttree i --minncells i --printdetail i\n");
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
  double tstep(0.01), tol(1e-3);
  double emass(0.511); //electron
  size_t npts(5000);
  bool draw(false), ttree(true);
  int minncells(10); // minimum # of hits
  //  double mine(90.0); // minimum energy to simulate
  // ttree variables
  TTree* cetree_;
  float targetde_, ipade_, trackerde_;
  float cee_, targete_, ipae_, trackere_;
  VEC3 cepos_, cemom_;
  float cet_;
  int npieces_, ntarget_, nipa_, ntrackerarcs_, ntrackercells_;
  int itrk_(0), kkstatus_, kkndof_, kknbf_, kknmat_, kknhit_, kkniter_;
  float kkchisq_, kkprob_;
  VEC3 kkentmom_, kkmidmom_, kkextmom_;
  // fit parameters
  KinKal::DVEC sigmas(0.5, 0.5, 0.5, 0.5, 0.002, 0.5); // expected parameter sigmas for loop helix
  double seedsmear(10.0);
  double dwt(1.0e6);
  unsigned maxniter(10);
//  Config::BFCorr bfcorr(Config::variable);
  Config::BFCorr bfcorr(Config::nocorr);
  Config::printLevel detail(Config::none);

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
    {"printdetail",     required_argument, 0, 'p' },
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
      case 'p' : detail = (Config::printLevel)atoi(optarg);
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
  // setup BField
  FileFinder filefinder;
  std::string fullfile = filefinder.fullFile(bfile);
  cout << "Building BField from file " << fullfile << endl;
//  AxialBFieldMap bfield(fullfile);
  UniformBFieldMap bfield(1.0);
//  GradientBFieldMap bfield(1.0,1.01, -1500,1500);

  bfield.print(cout);
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
  TH1F* ipade = new TH1F("ipde","IPA dE",100,0.0,0.5);
  TH1F* tarde = new TH1F("tarde","Target dE;dE (MeV)",100,0.001,5.0);
  TH1F* trkde = new TH1F("trkde","Tracker dE;dE (MeV)",100,0.001,2.0);
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
    cetree_->Branch("npieces",&npieces_,"npieces/I");
    cetree_->Branch("nipa",&nipa_,"nipa/I");
    cetree_->Branch("ipade",&ipade_,"ipade/F");
    cetree_->Branch("ntrackerarcs",&ntrackerarcs_,"ntrackerarcs/I");
    cetree_->Branch("ntrackercells",&ntrackercells_,"ntrackercells/I");
    cetree_->Branch("trackere",&trackere_,"trackere/F");
    cetree_->Branch("trackerde",&trackerde_,"trackerde/F");
    cetree_->Branch("targete",&targete_,"targete/F");
    cetree_->Branch("ipae",&ipae_,"ipae/F");
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
    targetde_ = ipade_ = trackerde_ = 0.0;
    targete_ = ipae_ = 0.0;
    nipa_ = ntrackerarcs_ = ntrackercells_ = -1;
    kkndof_ = kknbf_ = kknmat_ = kknhit_ = kkniter_ = -1;
    kkchisq_ = kkprob_ = -1.0;
    kkentmom_ = kkmidmom_ = kkextmom_ = VEC3();

    // generate a random CeEndpoint momentum; for now, use the endpoint energy, could use spectrum TODO
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
    double mom = sqrt(cee_*cee_ - emass*emass);
    cemom_ = VEC3(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost);
    ParticleState cestate(cepos_,cemom_,cet_,emass,-1);
    TimeRange range(cestate.time(),cestate.time()+1000.0); // replace hard-coded limit with a physical estimate FIXME
    auto bstart = bfield.fieldVect(cestate.position3());
    KTRAJ lhelix(cestate,bstart,range);
    //    cout << "Initial trajectory " << lhelix << endl;
    // initialize piecetraj
    PKTRAJ pktraj(lhelix);
    TimeRanges targetinters, ipainters, trackerinters;
    // extend through the target
    if(target.extendTrajectory(bfield,pktraj,targetinters)){
//      cout << "Extended to target " << targetinters.size() << endl;
      targete_ = pktraj.energy(pktraj.range().end());
      targetde_ = targete_ - cee_;
      ntarget_ = targetinters.size();
      // extend through the IPA
      if(ipa.extendTrajectory(bfield,pktraj,ipainters)){
//        cout << "Extended to ipa " << ipainters.size() << endl;
        nipa_ = ipainters.size();
        ipae_ = pktraj.energy(pktraj.range().end());
        ipade_ = ipae_ - targete_;
        // extend through tracker, and create the material xings and hits
        std::vector<double> htimes;
        std::vector<std::shared_ptr<Hit<KTRAJ>>> hits;
        std::vector<std::shared_ptr<ElementXing<KTRAJ>>> xings;
        double trackerpath(0.0);
        double speed = pktraj.speed(pktraj.range().end());
        tracker.simulateHits(bfield,pktraj,hits,xings,trackerinters,htimes);
        ntrackercells_ = htimes.size();
        for(auto const& inter : trackerinters) { trackerpath += speed*inter.range(); }
        if(ntrackercells_ > minncells){
//          cout << "Extended to tracker " << trackerinters.size() << endl;
          ntrackerarcs_ = trackerinters.size();
          double ke = cestate.energy() - cestate.mass();
          trackerde_ = -100*trackerEStar.dEIonization(ke)*tracker.density()*trackerpath; // unit conversion
          trackere_ = pktraj.energy(pktraj.range().end());
          tarde->Fill(targetde_);
          ipade->Fill(ipade_);
          trkde->Fill(trackerde_);
          trknc->Fill(ntrackercells_);
          // reconstruct track from hits and material info
          auto seedtraj = pktraj.nearestPiece(0.5*(htimes.front()+htimes.back()));
          for(size_t ipar=0; ipar < NParams(); ipar++){
            double perr = sigmas[ipar]*seedsmear;
            seedtraj.params().covariance()[ipar][ipar] = perr*perr;
//            seedtraj.params().parameters()[ipar] += tr_.Gaus(0.0,perr);
          }
          if(config.plevel_ > Config::none) {
            cout << "True traj " << pktraj.nearestPiece(trackerinters.front().begin()) << endl;
            cout << "Seed traj " << seedtraj << endl;
          }
          KKTRK kktrk(config,bfield,seedtraj,hits,xings);
          // fill fit information
          auto const& fstat = kktrk.fitStatus();
          kkstatus_ = fstat.status_;
          kkchisq_ = fstat.chisq_.chisq();
          kkprob_ = fstat.chisq_.probability();
          kkndof_ = fstat.chisq_.nDOF();
          kknbf_ = 0; kknhit_ = 0; kknmat_ = 0;
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
            auto const& fptraj = kktrk.fitTraj();
            auto entstate = fptraj.stateEstimate(fptraj.ztime(tracker.zMin()));
            auto midstate = fptraj.stateEstimate(fptraj.ztime(0.0));
            auto extstate = fptraj.stateEstimate(fptraj.ztime(tracker.zMax()));
            kkentmom_ = entstate.momentum3();
            kkmidmom_ = midstate.momentum3();
            kkextmom_ = extstate.momentum3();
          }
          cetree_->Fill();
          //
          if(draw){
            plhel.push_back(new TPolyLine3D(npts));
            plhel.back()->SetLineColor(icolor++%10);
            double tstart = pktraj.range().begin();
            double ts = pktraj.range().range()/(npts-1);
            VEC3 ppos;
            for(unsigned ipt=0;ipt<npts;ipt++){
              double t = tstart + ipt*ts;
              ppos = pktraj.position3(t);
              plhel.back()->SetPoint(ipt,ppos.X(),ppos.Y(),ppos.Z());
            }
          }
        } // particle hits the tracker
      } // particle stops in IPA
    } // particle exits the target going downstream
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
    rulers->SetAxisRange(bfield.zMin(),tracker.cylinder().zmax(),"Z");
    rulers->Draw();
    cetcan->Write();
  }

  cetracksfile.Write();
  cetracksfile.Close();

  // now draw
  mustopsfile->Close();
  return 0;
}

