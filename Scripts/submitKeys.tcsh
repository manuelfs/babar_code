#!/bin/tcsh

cd ..
srtpath `cat .current` $BFARCH
cd -

set q = long
set x = xlong
set d = `date +%y-%m-%d`
set logdir = AWG82/logfiles/results/keys/TimesFitPdf/$d/
mkdir -p $logdir

set sample = RAll
set Times = 1

bsub -q $q -o $logdir/logFit_01.log FitPdfSampleKEYS  1  80 240   1000 300 100 50  0.4   40 250  1.3   noRot  $sample  120 100 2     $Times # Sample 01
bsub -q $q -o $logdir/logFit_02.log FitPdfSampleKEYS  2  100 200  1000 300 100 40  0.4   60 180  1.5   noRot  $sample  110 100 2     $Times # Sample 02
bsub -q $q -o $logdir/logFit_03.log FitPdfSampleKEYS  3  80 200   1000 300 100 50  0.4   80 50   -1.5  noRot  $sample  120 100 2     $Times # Sample 03
bsub -q $q -o $logdir/logFit_04.log FitPdfSampleKEYS  4  120 200  1000 300 110 40  0.4   60 200  1.3   noRot  $sample  140 100 4     $Times # Sample 04
bsub -q $q -o $logdir/logFit_05.log FitPdfSampleKEYS  5  90 160   1000 300 80 50   0.4   50 200  1.5   noRot  $sample  120 100 3     $Times # Sample 05
bsub -q $q -o $logdir/logFit_06.log FitPdfSampleKEYS  6  80 160   1000 300 60 100  0.4   60 180  1.3   noRot  $sample  100 100 -5    $Times # Sample 06
bsub -q $q -o $logdir/logFit_07.log FitPdfSampleKEYS  7  80 160   1000 300 100 35  0.4   40 240  1.4   noRot  $sample  100 80 5      $Times # Sample 07
bsub -q $q -o $logdir/logFit_08.log FitPdfSampleKEYS  8  130 100  1000 300 120 80  0.4   70 80   1.5   noRot  $sample  40 150 5      $Times # Sample 08
bsub -q $x -o $logdir/logFit_09.log FitPdfSampleKEYS  9  30 500   1000 300 30 100  -0.9  30 100  2.0   noRot  $sample  230  35 1.1   $Times # Sample 09
bsub -q $x -o $logdir/logFit_10.log FitPdfSampleKEYS 10  30 400   1000 300 25 100  -0.9  25 100  2.0   noRot  $sample  220  35 1.6   $Times # Sample 10
bsub -q $x -o $logdir/logFit_11.log FitPdfSampleKEYS 11  30 520   1000 300 30 150  -0.9  30 150  2.0   noRot  $sample  230  35 1.1   $Times # Sample 11
bsub -q $x -o $logdir/logFit_12.log FitPdfSampleKEYS 12  35 470   1000 300 30 150  -0.9  30 160  2.0   noRot  $sample  260  40 1.6   $Times # Sample 12
bsub -q $x -o $logdir/logFit_13.log FitPdfSampleKEYS 13  70 280   1000 300 35 130  -0.4  45 130  1.9   noRot  $sample  250  40 3     $Times # Sample 13
bsub -q $q -o $logdir/logFit_14.log FitPdfSampleKEYS 14  70 230   1000 300 90 120  0.8   50 160  1.9   noRot  $sample  200  50 1     $Times # Sample 14
bsub -q $x -o $logdir/logFit_15.log FitPdfSampleKEYS 15  80 300   1000 300 70 100  -0.4  70 100  1.9   noRot  $sample  200  50 2.5   $Times # Sample 15
bsub -q $q -o $logdir/logFit_16.log FitPdfSampleKEYS 16  50 200   1000 300 150 100 1.0   25 350  1.5   noRot  $sample  100 100 -5    $Times # Sample 16
bsub -q $q -o $logdir/logFit_17.log FitPdfSampleKEYS 17  110 160  1000 300 70 100  0.4   60 110  1.7   noRot  $sample  120 140 5     $Times # Sample 17
bsub -q $q -o $logdir/logFit_18.log FitPdfSampleKEYS 18  100 160  1000 300 70 120  0.4   60 120  1.7   noRot  $sample  120 160 5     $Times # Sample 18
bsub -q $q -o $logdir/logFit_19.log FitPdfSampleKEYS 19  110 160  1000 300 70 100  0.4   50 110  1.7   noRot  $sample  120 140 5     $Times # Sample 19
bsub -q $q -o $logdir/logFit_20.log FitPdfSampleKEYS 20  110 160  1000 300 80 120  0.4   60 120  1.7   noRot  $sample  120 160 5     $Times # Sample 20

bsub -q $x -o $logdir/logFit_21.log FitPdfSampleKEYS 21  13 800   1000 300 100 30  0.4   13 150  1.9   noRot  $sample  120 100 0.4   $Times # Sample 21
bsub -q $x -o $logdir/logFit_22.log FitPdfSampleKEYS 22  13 800   1000 300 100 30  0.4   13 150  1.9   noRot  $sample  120 100 0.4   $Times # Sample 22
bsub -q $x -o $logdir/logFit_23.log FitPdfSampleKEYS 23  25 650   1000 300 100 40  0.4   25 150  1.9   noRot  $sample  120 100 0.4   $Times # Sample 23
bsub -q $x -o $logdir/logFit_24.log FitPdfSampleKEYS 24  16 800   1000 300 100 30  0.4   17 150  1.9   noRot  $sample  120 100 0.4   $Times # Sample 24
bsub -q $q -o $logdir/logFit_25.log FitPdfSampleKEYS 25  80 150   1000 300 100 70  1.4   20 280  -1.4  noRot  $sample   25 400 -0.2  $Times # Sample 25
bsub -q $q -o $logdir/logFit_26.log FitPdfSampleKEYS 26  80 150   1000 300 100 70  1.4   20 280  -1.4  noRot  $sample   25 450 -0.2  $Times # Sample 26
bsub -q $q -o $logdir/logFit_27.log FitPdfSampleKEYS 27  80 150   1000 300 100 70  1.4   20 280  -1.4  noRot  $sample   25 400 -0.2  $Times # Sample 27
bsub -q $q -o $logdir/logFit_28.log FitPdfSampleKEYS 28  80 160   1000 300 110 90  1.4   20 280  -1.4  noRot  $sample   25 450 -0.2  $Times # Sample 28
bsub -q $q -o $logdir/logFit_29.log FitPdfSampleKEYS 29  60 200   1000 300 100 30  0.3   50 50   2.0   noRot  $sample  150  50 1.6   $Times # Sample 29
bsub -q $q -o $logdir/logFit_30.log FitPdfSampleKEYS 30  60 180   1000 300 100 30  0.3   50 50   2.0   noRot  $sample  150  40 1.6   $Times # Sample 30
bsub -q $q -o $logdir/logFit_31.log FitPdfSampleKEYS 31  70 170   1000 300 100 30  0.4   40 70   1.9   noRot  $sample  150  40 1.6   $Times # Sample 31
bsub -q $q -o $logdir/logFit_32.log FitPdfSampleKEYS 32  70 300   1000 300 100 40  0.3   60 70   2     noRot  $sample  100 100 -5    $Times # Sample 32

bsub -q $q -o $logdir/logFit_41.log FitPdfSampleKEYS 41  40 350   1000 300 80 70   0.5   20 170  2.0   noRot  $sample  120  70 2.3   $Times # Sample 41
bsub -q $q -o $logdir/logFit_42.log FitPdfSampleKEYS 42  35 300   1000 300 100 40  0.5   30 130  1.9   noRot  $sample  180  35 2.3   $Times # Sample 42
bsub -q $q -o $logdir/logFit_43.log FitPdfSampleKEYS 43  110 300  1000 300 120 50  0.5   30 150  1.9   noRot  $sample  180  35 2.1   $Times # Sample 43
bsub -q $q -o $logdir/logFit_44.log FitPdfSampleKEYS 44  45 600   1000 300 120 70  1.4   25 240  2.0   noRot  $sample  180  35 2.1   $Times # Sample 44
bsub -q $q -o $logdir/logFit_45.log FitPdfSampleKEYS 45  90 180   1000 300 100 40  0.4   45 160  1.5   noRot  $sample  140  70 1.6   $Times # Sample 45
bsub -q $q -o $logdir/logFit_46.log FitPdfSampleKEYS 46  90 220   1000 300 100 30  0.4   45 180  1.5   noRot  $sample  140  90 1.6   $Times # Sample 46
bsub -q $q -o $logdir/logFit_47.log FitPdfSampleKEYS 47  90 160   1000 300 100 40  0.4   45 150  1.5   noRot  $sample  140  60 1.6   $Times # Sample 47
bsub -q $q -o $logdir/logFit_48.log FitPdfSampleKEYS 48  90 160   1000 300 100 40  0.4   45 150  1.5   noRot  $sample  140  60 1.6   $Times # Sample 48
bsub -q $q -o $logdir/logFit_51.log FitPdfSampleKEYS 51  80 200   1000 300 100 50  0.4   45 200  1.4   noRot  $sample  130  90 1.4   $Times # Sample 51
bsub -q $q -o $logdir/logFit_52.log FitPdfSampleKEYS 52  80 200   1000 300 100 35  0.4   45 200  1.4   noRot  $sample  140  90 1.4   $Times # Sample 52
bsub -q $q -o $logdir/logFit_53.log FitPdfSampleKEYS 53  90 200   1000 300 100 50  0.4   75 200  1.4   noRot  $sample  150  90 1.4   $Times # Sample 53
bsub -q $q -o $logdir/logFit_54.log FitPdfSampleKEYS 54  80 200   1000 300 100 50  0.4   65 200  1.4   noRot  $sample  160  90 1.4   $Times # Sample 54
bsub -q $x -o $logdir/logFit_55.log FitPdfSampleKEYS 55  110 200  1000 300 100 40  0.4   45 150  1.3   noRot  $sample  130  110 3    $Times # Sample 55
bsub -q $x -o $logdir/logFit_56.log FitPdfSampleKEYS 56  90 210   1000 300 100 30  0.3   45 180  1.3   noRot  $sample  140  140 3    $Times # Sample 56
bsub -q $x -o $logdir/logFit_57.log FitPdfSampleKEYS 57  90 210   1000 300 100 40  0.4   55 180  1.3   noRot  $sample  130  150 3    $Times # Sample 57
bsub -q $q -o $logdir/logFit_58.log FitPdfSampleKEYS 58  130 180  1000 300 120 40  0.4   45 100  1.3   noRot  $sample  130  140 4    $Times # Sample 58
bsub -q $q -o $logdir/logFit_61.log FitPdfSampleKEYS 61  110 150  1000 300 130 30  0.4   70 200  1.4   noRot  $sample  160  70 2     $Times # Sample 61
bsub -q $q -o $logdir/logFit_62.log FitPdfSampleKEYS 62  70 200   1000 300 100 30  0.3   50 220  1.4   noRot  $sample  120  90 2     $Times # Sample 62
bsub -q $q -o $logdir/logFit_63.log FitPdfSampleKEYS 63  100 180  1000 300 160 30  0.3   50 200  1.2   noRot  $sample  150  90 2     $Times # Sample 63
bsub -q $q -o $logdir/logFit_64.log FitPdfSampleKEYS 64  120 80   1000 300 45 130  -0.4  45 130  -1.8  noRot  $sample  100 100 -5    $Times # Sample 64
bsub -q $q -o $logdir/logFit_65.log FitPdfSampleKEYS 65  70 130   1000 300 100 40  0.35  50 200  1.2   noRot  $sample  130  80 2     $Times # Sample 65
bsub -q $q -o $logdir/logFit_66.log FitPdfSampleKEYS 66  70 170   1000 300 100 40  0.3   40 230  1.3   noRot  $sample  150  100 2    $Times # Sample 66
bsub -q $q -o $logdir/logFit_67.log FitPdfSampleKEYS 67  70 150   1000 300 100 40  0.3   40 200  1.3   noRot  $sample  150  100 2    $Times # Sample 67
bsub -q $q -o $logdir/logFit_68.log FitPdfSampleKEYS 68  110 100  1000 300 80 50   0.3   80 120  1.0   noRot  $sample  100 100 -5    $Times # Sample 68
							   
