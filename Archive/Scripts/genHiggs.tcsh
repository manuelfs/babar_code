#!/bin/tcsh

if ( $# < 1 ) then
    echo "The format is: sendHiggs.tcsh sample"
    printf "\n\n"
    exit -1
endif

set nSample = $1
set Sample  = $nSample
if ( $nSample < 100 ) then
    set Sample  = 0$Sample
endif
if ( $nSample < 10 ) then
    set Sample  = 0$Sample
endif

GenPdfSamples RAllHigx${Sample} no_YES

set q = medium
set x = xlong
set d = `date +%y-%m-%d`
set logdir = AWG82/logfiles/results/keys/FitPdf/${d}_${Sample}/
mkdir -p $logdir

set sample = RAllHigx${Sample}
set Times = 1

bsub -q $q -o $logdir/logFit_01.log FitPdfSampleKEYS  1  91      175	  1000 300 100 50  0.4   50 250  1.3   noRot  $sample  130 100 2     $Times # Sample 01
bsub -q $q -o $logdir/logFit_02.log FitPdfSampleKEYS  2  91      102	  1000 300 100 40  0.4   60 180  1.5   noRot  $sample  110 100 2     $Times # Sample 02
bsub -q $q -o $logdir/logFit_03.log FitPdfSampleKEYS  3  123     114	  1000 300 100 50  0.4   80 50   -1.5  noRot  $sample  120 100 2     $Times # Sample 03
bsub -q $q -o $logdir/logFit_04.log FitPdfSampleKEYS  4  83      151	  1000 300 110 40  0.4   60 200  1.3   noRot  $sample  175 100 4     $Times # Sample 04
bsub -q $q -o $logdir/logFit_05.log FitPdfSampleKEYS  5  108     97	  1000 300 80 50   0.4   50 200  1.5   noRot  $sample  120 100 3     $Times # Sample 05
bsub -q $q -o $logdir/logFit_06.log FitPdfSampleKEYS  6  167     101	  1000 300 60 100  0.4   60 180  1.3   noRot  $sample  100 100 -5    $Times # Sample 06
bsub -q $q -o $logdir/logFit_07.log FitPdfSampleKEYS  7  115     96	  1000 300 100 35  0.4   40 240  1.4   noRot  $sample  100 80 5      $Times # Sample 07
bsub -q $q -o $logdir/logFit_08.log FitPdfSampleKEYS  8  240     58	  1000 300 120 80  0.4   70 80   1.5   noRot  $sample  40 150 5      $Times # Sample 08

printf "\n\n"
echo babar_code/Scripts/mvHiggs.tcsh ${Sample}
@ newSample = $nSample + 5
echo babar_code/Scripts/genHiggs.tcsh ${newSample}
printf "\n\n"
