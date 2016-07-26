#!/bin/tcsh


if ( $# < 1 ) then
    set Tag = RAll
else
    set Tag = $1
endif

set folder = ~/releases/Ntuple51/workdir/keys/root/fitmuNewx100/
foreach file ( $folder/*.root )
    KeysPlot $file ${Tag} curve
end


