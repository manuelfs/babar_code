#!/bin/tcsh

# Runs Mazur's gscalc

if ( $# < 0 ) then
    printf "\n\n"
    echo "The format is: sendgscalc.tcsh nJobs"
    printf "\n\n"
    exit -1
endif
if ( $# > 0 ) then
    @ nJobs = $1
else
    @ nJobs = 10
endif

set logdir = babar_code/NewPhysics/logfiles
rm $logdir/*
@ index = 0
@ delta = ( 2000 + $nJobs - 1 ) / $nJobs

while ( $index < 2000 )
    @ Nend = $index + $delta
    if ( $Nend > 2000 ) @ Nend = 2000
    bsub -q medium -o $logdir/loggscalc_${index}_${Nend}.log babar_code/NewPhysics/gscalc.tcsh $index $Nend
    @ index += $delta
end
