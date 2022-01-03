//
// struct for summarizing track fit information
//
namespace TrackToy {
  struct TrkInfo {
    unsigned ncells, narcs, ntrkhits, ncalohits; // hit counting
    int kkstatus, kkndof, kknbf, kknmat, kknhit, kkniter;
    float kkchisq, kkprob;
    void reset() {
      ncells = narcs = ntrkhits = ncalohits = 0;
      kkstatus = kkndof = kknbf = kknmat = kknhit = kkniter = -1;
      kkchisq = kkprob = -1.0;
    }
  };
}
