#!/bin/tcsh

if ( $# < 1 ) then
    echo "The format is: sendHiggs.tcsh sample"
    printf "\n\n"
    exit -1
endif

set digit = $1
#set Samples = $1
#if ( $AllSamples == a0 ) set Samples = ( 010 020 030 040 050 060 070 080 090 100 )
#if ( $AllSamples == a5 ) set Samples = ( 005 015 025 035 045 055 065 075 085 095 )
set Samples = ( 00 01 02 03 04 05 06 07 08 09 10 )

foreach sample ( $Samples )
    #mkdir keys/root/fitHigx${sample}
    #cpf keys/root/fitNewx100/* keys/root/fitHigx${sample}
    MultiEventMixer AWG82/ntuples/small/RAll_RunAll.root Hig yes_YES_mES ${sample}${digit}
end
