#!/bin/tcsh

if ( $# < 1 ) then
    echo "The format is: sendHiggs.tcsh sample"
    printf "\n\n"
    exit -1
endif

set digit = $1
set Samples = ( 00 01 02 03 04 05 06 07 08 09 10 )

#if ( $digit == 0 ) then
#    cpf AWG82/ntuples/small/FitRAllNewx100_RunAll.root AWG82/ntuples/small/FitRAllHigx000_RunAll.root
#endif
foreach sample ( $Samples )
    set tBmH = ${sample}${digit}
    #if ( $tBmH == 000 ) then
    #	set tBmH = 100
    #endif
    #mkdir keys/root/fitHigx${tBmH}
    #cpf keys/root/fitNewx100/* keys/root/fitHigx${tBmH}
    set logfile = AWG82/logfiles/Higgs/logMultiHiggs_${tBmH}.log
    rm $logfile
    bsub -q medium -o $logfile MultiEventMixer AWG82/ntuples/small/RAll_RunAll.root Hig yes_YES_mES ${tBmH}
end
