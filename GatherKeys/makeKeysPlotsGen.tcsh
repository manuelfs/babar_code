#!/bin/tcsh


if ( $# < 1 ) then
    set Tag = All
else
    set Tag = $1
endif

set folder = ~/releases/Ntuple51/workdir/keys/root/fitSep/
foreach file ( $folder/*.root )
    KeysPlot $file ${Tag} diff 25 -2 10
end
set folder = ~/releases/Ntuple51/workdir/keys/root/fitComb/
foreach file ( $folder/*.root )
    KeysPlot $file ${Tag} diff 25 -3 10
end


