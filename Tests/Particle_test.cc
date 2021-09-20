//
//  Test program for particle propagation in a BField
//
#include "KinKal/General/AxialBFieldMap.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "TrackToy/General/FileFinder.hh"
#include "TFile.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
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

using namespace TrackToy;
using namespace std;

void print_usage() {
  printf("Usage: ParticleTest --pfile s --bfile s --zmax f --tol f --minmass f --npts i --ntrks i\n");
}

int main(int argc, char **argv) {
  using KTRAJ=KinKal::LoopHelix;
  using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
  size_t npts(1000);
  int ntrks(-1);
  string pfile, bfile;
  double zmax(-3500), tol(1e-3);
  double minmass(100.0);

  static struct option long_options[] = {
    {"pfile",     required_argument, 0, 'f' },
    {"bfile",     required_argument, 0, 'F' },
    {"tol",     required_argument, 0, 't' },
    {"minmass",     required_argument, 0, 'm' },
    {"zmax",     required_argument, 0, 'z'  },
    {"npts",     required_argument, 0, 'N'  },
    {"ntrks",     required_argument, 0, 'n'  },
    {NULL, 0,0,0}
  };
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'f' : pfile = string(optarg);
                 break;
      case 'F' : bfile = string(optarg);
                 break;
      case 't' : tol = atof(optarg);
                 break;
      case 'm' : minmass = atof(optarg);
                 break;
      case 'z' : zmax = atof(optarg);
                 break;
      case 'N' : npts = atoi(optarg);
                 break;
      case 'n' : ntrks = atoi(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  // not sure why this is necessary...
  gSystem->Load("lib/libTests.dylib");
  gSystem->Load("libGeom");
  if(pfile.size()==0){
    cout << "No input pfile specified: terminating" << endl;
    return 1;
  }
  // open the input particle pfile
  TFile* ppfile = TFile::Open(pfile.c_str(),"READ");
  // find the TTree in the pfile
  TDirectory* td = (TDirectory*)ppfile->Get("StepPointMCDumper");
  TTreeReader reader("nt",td);
  TTreeReaderValue<KinKal::ParticleState> pstate(reader, "particle");
  TTree* ptree = (TTree*)td->Get("nt");
  cout << "Particle TTree has " << ptree->GetEntries() << " Entries" << endl;
  // setup BField
  FileFinder filefinder;
  std::string fullfile = filefinder.fullFile(bfile);
  KinKal::AxialBFieldMap axfield(fullfile);
  cout << "axial field between " << axfield.zMin() << " and " << axfield.zMax() << " with " << axfield.field().size()
    << " field values from "  << axfield.field().front() << " to "  << axfield.field().back() << endl;
  TFile ptestpfile("ParticleTest.root","RECREATE");
  std::vector<TPolyLine3D*> plhel;
  int icolor(kBlue);
  zmax = std::min(zmax,axfield.zMax());
  int itrk(0);
  while (reader.Next() && (ntrks < 0 || itrk < ntrks)) {
    // select by mass
    if(pstate->mass() > minmass){
    itrk++;
//    cout << "pstate status " <<  pstate.GetSetupStatus() << endl;
//    cout << "Read particle with position " << pstate->position3() << " time " << pstate->time() << " momentum " << pstate->momentum3() << endl;
//    cout << "pVz = " << pstate->velocity().Z() << " tmax " << (zmax-pstate->position3().Z())/pstate->velocity().Z() << endl;

    KinKal::TimeRange range(pstate->time(),pstate->time()+(zmax-pstate->position3().Z())/pstate->velocity().Z());
//    cout << " initial range " << range << endl;
    // create a loop helix from this
    auto bstart = axfield.fieldVect(pstate->position3());
    KTRAJ lhelix(*pstate,bstart,range);
//    cout << "Initial trajectory " << lhelix << endl;
    // initialize piecetraj
    PKTRAJ phelix(lhelix);
    auto pos = pstate->position3();
    // extend to the end of the BField range
    while(pos.Z() < zmax){
      range.begin() = axfield.rangeInTolerance(phelix.back(),range.begin(),tol);
      if(range.begin() < range.end()){
        // Predict new position and momentum at this end, making linear correction for BField effects
        auto endstate = phelix.back().state(range.begin());
        pos = endstate.position3();
        auto bend = axfield.fieldVect(pos);
        KTRAJ endhelix(endstate,bend,range);
        phelix.append(endhelix);
        //      cout << "appended helix at point " << pos << " time " << range.begin() << endl;
      } else
        break;
    }
    double largest, average;
    size_t igap;
    phelix.gaps(largest, igap, average);
//    cout << "Particle piece traj with " << phelix.pieces().size() << " pieces and largest gap = "
//      << largest << " average gap = " << average << endl;
    plhel.push_back(new TPolyLine3D(npts));
    plhel.back()->SetLineColor(icolor++%10);
    double tstart = phelix.range().begin();
    double ts = phelix.range().range()/(npts-1);
    KinKal::VEC3 ppos;
    for(unsigned ipt=0;ipt<npts;ipt++){
      double t = tstart + ipt*ts;
      ppos = phelix.position3(t);
      plhel.back()->SetPoint(ipt,ppos.X(),ppos.Y(),ppos.Z());
    }
    }
  }
// now draw
  TCanvas* ptcan = new TCanvas("Particle");
// Draw target

  TTUBE* target= new TTUBE("target","target","void",21.5,75.0,800.0);
  target->SetLineColor(kBlack);
  target->SetLineWidth(4);
  target->SetFillColorAlpha(kBlack, 0.5);
  target->Draw();
  TTUBE *DS  = new TTUBE("DS","DS","void",0.0,800,4000);
  DS->SetVisibility(0);
  DS->Draw();
  TNode *node1 = new TNode("NODE1","NODE1","DS");
  node1->cd();
  TNode* tnode = new TNode("targetnode","targetnode","target",0.0,0.0,-4300);
  tnode->Draw();
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
  rulers->SetAxisRange(axfield.zMin(),zmax,"Z");
  rulers->Draw();

  ptcan->Write();
  ptestpfile.Write();
  ptestpfile.Close();
  ppfile->Close();
  return 0;
}

