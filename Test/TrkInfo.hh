//
// struct for summarizing track fit information
//
namespace TrackToy {
  struct TrkInfo {
    unsigned ncells, narcs, ntrkhits, ncalohits; // hit counting
    int kkstatus, kkndof, kknbf, kknmat, kknhit, kkniter;
    int kknactive, kknnull;
    float kkchisq, kkprob;
    void reset() {
      ncells = narcs = ntrkhits = ncalohits = 0;
      kkstatus = -1;
      kkndof = kknbf = kknmat = kknhit = kkniter = 0;
      kknactive = kknnull = 0;
      kkchisq = kkprob = -1.0;
    }
  };
}
