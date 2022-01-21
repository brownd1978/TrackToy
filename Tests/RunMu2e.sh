#
#  Only need to run MuStops once
bin/MuStops --mubeam MDC2020n_10pc --targetfile Data/Mu2eTarget.dat >& MDC2020N_10pc_mustops.log
nohup time bin/Tracks --process CeMinus --rmue 1e-16 --ntrks 100000 --mustopsfile MDC2020n_10pc --npot 3.6e20 --tfile Mu2eMCAmbig --trackerfile Data/Mu2eTracker.dat --fitschedule Data/Schedule.txt >& Mu2eMCAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --mustopsfile MDC2020n_10pc --npot 3.6e20 --tfile Mu2eMCAmbig --trackerfile Data/Mu2eTracker.dat --fitschedule Data/Schedule.txt >& Mu2eMCAmbigDIO.log &
nohup time bin/Tracks --process CeMinus --rmue 1e-16 --ntrks 100000 --mustopsfile MDC2020n_10pc --npot 3.6e20 --tfile Mu2eDriftAmbig --trackerfile Data/Mu2eTracker.dat --fitschedule Data/Schedule_driftfit.txt >& Mu2eDriftAmbigCeM.log &
nohup time bin/Tracks --process FlatDIO --ntrks 100000 --mustopsfile MDC2020n_10pc --npot 3.6e20 --tfile Mu2eDriftAmbig --trackerfile Data/Mu2eTracker.dat --fitschedule Data/Schedule_driftfit.txt >& Mu2eDriftAmbigDIO.log &

