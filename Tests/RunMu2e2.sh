nohup time bin/Tracks --process CeMinus --ntrks 100000 --tfile Mu2e2MCAmbig --trackerfile Data/Mu2e2Tracker.dat >& Mu2e2MCAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --tfile Mu2e2_MCAmbig --trackerfile Data/Mu2e2Tracker.dat >& Mu2e2_MCAmbigDIO.log &
nohup time bin/Tracks --process CeMinus --ntrks 100000 --tfile Mu2e2DriftAmbig --trackerfile Data/Mu2e2Tracker.dat --fitschedule Data/Schedule_driftfit.txt >& Mu2e2DriftAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --tfile Mu2e2DriftAmbig --trackerfile Data/Mu2e2Tracker.dat --fitschedule Data/Schedule_driftfit.txt >& Mu2e2DriftAmbigDIO.log &

