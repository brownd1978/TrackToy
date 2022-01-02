//
// struct for summarizing track fit information
//
namespace TrackToy {
  struct TrkInfo {
    float wt_; // weight, coming from generator and general parameters
    unsigned ncells_, narcs_, ntrkhits_, ncalohits_; // hit counting
    int kkstatus_, kkndof_, kknbf_, kknmat_, kknhit_, kkniter_;
    float kkchisq_, kkprob_;
    void reset() {
      ncells_ = narcs_ = ntrkhits_ = ncalohits_ = 0;
      kkstatus_ = kkndof_ = kknbf_ = kknmat_ = kknhit_ = kkniter_ = -1;
      kkchisq_ = kkprob_ = -1.0;
    }
  };
}
