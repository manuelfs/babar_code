#!/bin/tcsh

if ( $# > 0 ) then
    set Multi = $1
else
    set Multi = 2
endif

set d = `date +%y-%m-%d`
set logdir = AWG82/logfiles/results/keys/Pulls/$d/
mkdir -p $logdir

@ nFiles = 2 * $Multi
@ iFile = 1
while ( $iFile <= $nFiles)
    set Tag = R24Rests_${iFile}_${nFiles}
    bsub -q medium -o $logdir/log${Tag}.log DonutFitFinal FitAll/configRAll.txt 0 1 $Tag $Tag
    @ iFile++
end
@ nFiles = 7 * $Multi
@ iFile = 1
while ( $iFile <= $nFiles)
    set Tag = R26Rests_${iFile}_${nFiles}
    bsub -q medium -o $logdir/log${Tag}.log DonutFitFinal FitAll/configRAll.txt 0 1 $Tag $Tag
    @ iFile++
end




