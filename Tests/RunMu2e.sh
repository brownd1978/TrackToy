nohup time bin/Tracks --process CeMinus --ntrks 100000 --tfile Mu2eMCAmbig --trackerfile Data/Mu2eTracker.dat >& Mu2eMCAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --tfile Mu2e_MCAmbig --trackerfile Data/Mu2eTracker.dat >& Mu2e_MCAmbigDIO.log &
nohup time bin/Tracks --process CeMinus --ntrks 100000 --tfile Mu2eDriftAmbig --trackerfile Data/Mu2eTracker.dat --fitschedule Data/Schedule_driftfit.txt >& Mu2eDriftAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --tfile Mu2eDriftAmbig --trackerfile Data/Mu2eTracker.dat --fitschedule Data/Schedule_driftfit.txt >& Mu2eDriftAmbigDIO.log &

