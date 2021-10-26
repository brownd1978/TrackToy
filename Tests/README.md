
# TrackToy/Tests
To run these test programs, you must have built both KinKal and TrackToy, and run:
> source setup.sh
in your build directory

## MuStops_test;
Find muon stopping positions given the target description and a file of the beam particles, produced by the StepPointMCDumper module in Mu2e/Offline
This produces the TFile and TTree used in subsequent tests that start from muon stopping positions
## CeTracks_test
generate Ce particles from stopped muon positions and propagate them through the passive materials and tracker
