#!/bin/tcsh

if ( $# > 0 ) then
    set Multi = $1
else
    set Multi = 1
endif

@ nFiles = 2 * $Multi
@ iFile = 1
while ( $iFile <= $nFiles)
    EventMixer AWG82/ntuples/Defsmall/R24Rests_${iFile}_${nFiles}_RunAll.root
    @ iFile++
end
@ nFiles = 7 * $Multi
@ iFile = 1
while ( $iFile <= $nFiles)
    EventMixer AWG82/ntuples/Defsmall/R26Rests_${iFile}_${nFiles}_RunAll.root
    @ iFile++
end




