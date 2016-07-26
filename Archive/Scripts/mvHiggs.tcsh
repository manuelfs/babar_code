#!/bin/tcsh

if ( $# < 1 ) then
    echo "The format is: sendHiggs.tcsh sample"
    printf "\n\n"
    exit -1
endif

set Sample = $1
set Indices = ( 1 2 3 4 5 6 7 8 )

foreach index ( $Indices )
    cpf keys/root/fitFinal/pdfKeys_${index}_Fit.root keys/root/fitHigx${Sample}/
end
