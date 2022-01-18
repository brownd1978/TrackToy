nohup time bin/Tracks --process CeMinus --ntrks 1000 --tfile TestFit --trackerfile Data/Mu2eTracker.dat --fitschedule Data/Schedule_driftfit.txt >& TestFitCeM.log &
nohup time bin/Tracks --process CeMinus --ntrks 1000 --tfile TestExt --trackerfile Data/Mu2eTracker.dat --bfit 0 --fitschedule Data/Schedule_seedfit.txt --extschedule Data/Schedule_driftextend.txt >& TestExtCeM.log &

