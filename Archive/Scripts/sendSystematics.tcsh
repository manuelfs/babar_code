#!/bin/tcsh

# Sends the variations for the systematics in AWG82/systematics/BF/

if ( $# < 1 ) then
    printf "\n\n"
    echo "The format is: sendgSystematics.tcsh nJobs nBunch"
    printf "\n\n"
    exit -1
endif

@ nJobs = $1
if ( $# > 1 ) then
    @ nBunch = $2
else
    @ nBunch = 10
endif

set logdir = AWG82/logfiles/systematics/BF/
rm $logdir/*
@ index = 0

while ( $index < $nJobs + 1 )
    @ Nend = $index + $nBunch - 1
    if ( $Nend > $nJobs ) @ Nend = $nJobs
    set logfile = $logdir/logSysBF_${index}_${Nend}.log
    if( ! -e $logfile ) bsub -q long -o $logfile Systematics $index $Nend
    @ index += $nBunch
end
