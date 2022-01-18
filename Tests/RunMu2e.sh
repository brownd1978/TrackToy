nohup time bin/Tracks --process CeMinus --ntrks 100000 --mustopsfile MDC2020n_10pc --tfile Mu2eMCAmbig --trackerfile Data/Mu2eTracker.dat --fitschedule Data/Schedule.txt >& Mu2eMCAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --mustopsfile MDC2020n_10pc --tfile Mu2eMCAmbig --trackerfile Data/Mu2eTracker.dat --fitschedule Data/Schedule.txt >& Mu2eMCAmbigDIO.log &
nohup time bin/Tracks --process CeMinus --ntrks 100000 --mustopsfile MDC2020n_10pc --tfile Mu2eDriftAmbig --trackerfile Data/Mu2eTracker.dat --fitschedule Data/Schedule_driftfit.txt >& Mu2eDriftAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --mustopsfile MDC2020n_10pc --tfile Mu2eDriftAmbig --trackerfile Data/Mu2eTracker.dat --fitschedule Data/Schedule_driftfit.txt >& Mu2eDriftAmbigDIO.log &

