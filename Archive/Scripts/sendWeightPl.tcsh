#!/bin/tcsh

cd ..
srtpath `cat .current` $BFARCH
cond24boot09
cd -
printf "\n"


set logdir = AWG82/logfiles/results/WeightPl
mkdir -p $logdir
rm ${logdir}/*log

set AllSample = (1 2 3 4)
foreach sample ( $AllSample )
    bsub -o ${logdir}/logWeightPl_${sample}.log  -q long  WeightPl 90 $sample 1
end

#set doText = t
#if( doText == t)
#    cat keys/weights/ContiWeightText_1_90.txt >! babar_code/Reweight/wContiPl
#    cat keys/weights/ContiWeightText_2_90.txt >> babar_code/Reweight/wContiPl
#    cat keys/weights/ContiWeightText_3_90.txt >> babar_code/Reweight/wContiPl
#    cat keys/weights/ContiWeightText_4_90.txt >> babar_code/Reweight/wContiPl
#endif

printf "\n"
