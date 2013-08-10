#!/bin/tcsh

# Runs Mazur's gscalc

if ( $# < 1 ) then
    printf "\n\n"
    echo "The format is: gscalc.tcsh iBegin [iEnd]"
    printf "\n\n"
    exit -1
endif
if ( $# > 0 ) then
    @ index = $1
else
    @ index = 1
endif
if ( $# > 1 ) then
    @ iEnd = $2
else
    @ iEnd = 1000
endif 

while ( $index < $iEnd )
    gscalc $index 1000
    @ index++
end
