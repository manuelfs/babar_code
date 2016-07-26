#!/bin/tcsh

cd ..
srtpath `cat .current` $BFARCH
cd -
printf "\n"

set logdir = AWG82/logfiles/results/OptiBDT
mkdir -p $logdir

set AllSample = (sig dss)
set AllType = (Dl Comb 2D)
set AllAge  = (new old)
set AllNtu  = (Add_R24Rest Add_R24Rest2)

set AllSample = (sepsig sepdss)
set AllType = (Dl Comb 2D)
set AllAge  = (new)
set AllNtu  = (Add_RAll)

foreach sample ( $AllSample )
    foreach type ( $AllType )
	foreach ntu ( $AllNtu )
	    foreach age ( $AllAge )
	    	set tag = old
	    	if ( $age == new ) set tag = ""
	    	set typeSample = ${sample}${type}
	    	set logfile = ${logdir}/${typeSample}_${ntu}.log
		set range = "0 0.85"
	    	if ( $sample == dss ) set range = "-0.7 0.1"
		set nbins = 85
	    	if ( $type == 2D ) set nbins = 35
		rm $logfile
	    	bsub -o $logfile -q long OptiBDT $typeSample $nbins $range ${tag}${ntu}  
	    end
	end
    end
end

printf "\n"
