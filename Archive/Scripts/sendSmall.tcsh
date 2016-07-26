#!/bin/tcsh

cd ..
srtpath `cat .current` $BFARCH
cond24boot09
cd -
printf "\n"

set logdir = AWG82/logfiles/ntuples/small
mkdir -p $logdir
rm ${logdir}/*log
set Tag = ""
if ( $# == 1 ) then
    set Tag = $1
endif


bsub -o ${logdir}/OffPeak_${Tag}.log           -q long SmallNtuple OffPeak_mES  All Whole      1				 $Tag  
bsub -o ${logdir}/Data_${Tag}.log              -q long SmallNtuple Data_mES     All Whole      1			         $Tag  
bsub -o ${logdir}/ccbar_${Tag}.log             -q long SmallNtuple ccbar_mES    All Whole      babar_code/Reweight/wTotal.txt  $Tag  
bsub -o ${logdir}/uds_${Tag}.log               -q long SmallNtuple uds_mES      All Whole      babar_code/Reweight/wTotal.txt  $Tag  
bsub -o ${logdir}/R24MVA${Tag}.log             -q long SmallNtuple R24_mES      All MVA        babar_code/Reweight/wTotal.txt  $Tag  
bsub -o ${logdir}/R24Rests_1_2${Tag}.log       -q long SmallNtuple R24_mES      All Rest  babar_code/Reweight/wTotal.txt  $Tag  
bsub -o ${logdir}/R24Rests_2_2${Tag}.log       -q long SmallNtuple R24_mES      All Rest2  babar_code/Reweight/wTotal.txt  $Tag  
bsub -o ${logdir}/R26Rests_1_7_${Tag}.log      -q long SmallNtuple R26_mES      All Rests_1_7  babar_code/Reweight/wTotal.txt  $Tag  
bsub -o ${logdir}/R26Rests_2_7_${Tag}.log      -q long SmallNtuple R26_mES      All Rests_2_7  babar_code/Reweight/wTotal.txt  $Tag  
bsub -o ${logdir}/R26Rests_3_7_${Tag}.log      -q long SmallNtuple R26_mES      All Rests_3_7  babar_code/Reweight/wTotal.txt  $Tag  
bsub -o ${logdir}/R26Rests_4_7_${Tag}.log      -q long SmallNtuple R26_mES      All Rests_4_7  babar_code/Reweight/wTotal.txt  $Tag  
bsub -o ${logdir}/R26Rests_5_7_${Tag}.log      -q long SmallNtuple R26_mES      All Rests_5_7  babar_code/Reweight/wTotal.txt  $Tag  
bsub -o ${logdir}/R26Rests_6_7_${Tag}.log      -q long SmallNtuple R26_mES      All Rests_6_7  babar_code/Reweight/wTotal.txt  $Tag  
bsub -o ${logdir}/R26Rests_7_7_${Tag}.log      -q long SmallNtuple R26_mES      All Rests_7_7  babar_code/Reweight/wTotal.txt  $Tag  

#bsub -o ${logdir}/R24Rests_1_4${Tag}.log       -q long SmallNtuple R24      All Rests_1_4  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R24Rests_2_4${Tag}.log       -q long SmallNtuple R24      All Rests_2_4  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R24Rests_3_4${Tag}.log       -q long SmallNtuple R24      All Rests_3_4  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R24Rests_4_4${Tag}.log       -q long SmallNtuple R24      All Rests_4_4  babar_code/Reweight/wTotal.txt  $Tag  
#
#bsub -o ${logdir}/R26Rests_1_14_${Tag}.log     -q long SmallNtuple R26      All Rests_1_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_2_14_${Tag}.log     -q long SmallNtuple R26      All Rests_2_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_3_14_${Tag}.log     -q long SmallNtuple R26      All Rests_3_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_4_14_${Tag}.log     -q long SmallNtuple R26      All Rests_4_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_5_14_${Tag}.log     -q long SmallNtuple R26      All Rests_5_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_6_14_${Tag}.log     -q long SmallNtuple R26      All Rests_6_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_7_14_${Tag}.log     -q long SmallNtuple R26      All Rests_7_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_8_14_${Tag}.log     -q long SmallNtuple R26      All Rests_8_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_9_14_${Tag}.log     -q long SmallNtuple R26      All Rests_9_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_10_14_${Tag}.log    -q long SmallNtuple R26      All Rests_10_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_11_14_${Tag}.log    -q long SmallNtuple R26      All Rests_11_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_12_14_${Tag}.log    -q long SmallNtuple R26      All Rests_12_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_13_14_${Tag}.log    -q long SmallNtuple R26      All Rests_13_14  babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/R26Rests_14_14_${Tag}.log    -q long SmallNtuple R26      All Rests_14_14  babar_code/Reweight/wTotal.txt  $Tag  
#
#bsub -o ${logdir}/R26_${Tag}.log               -q long SmallNtuple R26      All Whole      babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/cocktail_${Tag}.log          -q long SmallNtuple cocktail All Whole      babar_code/Reweight/wTotal.txt  $Tag  
#bsub -o ${logdir}/Dpipi_${Tag}.log             -q long SmallNtuple Dpipi    All Whole      babar_code/Reweight/wTotal.txt  $Tag  

printf "\n"
