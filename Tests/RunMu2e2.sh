bin/MuStops --mubeam Mu2eIIa2 --targetfile Data/Mu2e2Target.dat >& Mu2eIIa2_mustops.log
nohup time bin/Tracks --process CeMinus --rmue 3e-17 --ntrks 100000 --mustopsfile Mu2eIIa2 --npot 1.3e22 --tfile Mu2e2MCAmbig --trackerfile Data/Mu2e2Tracker.dat --fitschedule Schedule_MCAmbig.txt >& Mu2e2MCAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --mustopsfile Mu2eIIa2 --npot 1.3e22 --tfile Mu2e2MCAmbig --trackerfile Data/Mu2e2Tracker.dat --fitschedule Schedule_MCAmbig.txt >& Mu2e2MCAmbigFlatDIO.log &
nohup time bin/Tracks --process CeMinus --rmue 3e-17 --ntrks 100000 --mustopsfile Mu2eIIa2 --npot 1.3e22 --tfile Mu2e2DriftAmbig --trackerfile Data/Mu2e2Tracker.dat --fitschedule Schedule_driftfit.txt >& Mu2e2DriftAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --mustopsfile Mu2eIIa2 --npot 1.3e22 --tfile Mu2e2DriftAmbig --trackerfile Data/Mu2e2Tracker.dat --fitschedule Schedule_driftfit.txt >& Mu2e2DriftAmbigFlatDIO.log &
nohup time bin/Tracks --process DIO --endrange 3.0 --ntrks 100000 --mustopsfile Mu2eIIa2 --npot 1.3e22 --tfile Mu2e2DriftAmbig --trackerfile Data/Mu2e2Tracker.dat --fitschedule Schedule_driftfit.txt >& Mu2e2DriftAmbigDIO.log &

