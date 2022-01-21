bin/MuStops --mubeam Mu2eIIa2 --targetfile Data/Mu2e2Target.dat >& Mu2eII2a_mustops.log
nohup time bin/Tracks --process CeMinus --rmue 1e-17 --ntrks 100000 --mustopsfile MDC2020n_10pc --npot 1.3e22 --tfile Mu2e2MCAmbig --trackerfile Data/Mu2e2Tracker.dat >& Mu2e2MCAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --mustopsfile MDC2020n_10pc --npot 1.3e22 --tfile Mu2e2MCAmbig --trackerfile Data/Mu2e2Tracker.dat >& Mu2e2MCAmbigDIO.log &
nohup time bin/Tracks --process CeMinus --rmue 1e-17 --ntrks 100000 --mustopsfile MDC2020n_10pc --npot 1.3e22 --tfile Mu2e2DriftAmbig --trackerfile Data/Mu2e2Tracker.dat --fitschedule Data/Schedule_driftfit.txt >& Mu2e2DriftAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --mustopsfile MDC2020n_10pc --npot 1.3e22 --tfile Mu2e2DriftAmbig --trackerfile Data/Mu2e2Tracker.dat --fitschedule Data/Schedule_driftfit.txt >& Mu2e2DriftAmbigDIO.log &

