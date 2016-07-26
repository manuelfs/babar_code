#!/bin/tcsh


set folder = ~/releases/mike50/workdir/keys/root/fitSep/

foreach file ( $folder/*.root )
    KeysPlot $file cocktailTau
end


