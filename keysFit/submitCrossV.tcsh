#!/usr/local/bin/tcsh

cd ..
srtpath `cat .current` $BFARCH
anal50boot
cd -

set q = xlong
set x = xxl
set d = `date +%y-%m-%d`
set logdir = AWG82/logfiles/results/keys/CrossV/$d/
mkdir -p $logdir

bsub -q $x -o $logdir/logCrossV_0.log CrossValidation 0  6 97 118 5 132 161     # Sample 0
bsub -q $x -o $logdir/logCrossV_1.log CrossValidation 1  9 75 92 6 124 151      # Sample 1
bsub -q $x -o $logdir/logCrossV_2.log CrossValidation 2  6 81 99 3 132 161      # Sample 2
bsub -q $x -o $logdir/logCrossV_3.log CrossValidation 3  4 68 82 3 135 165     # Sample 3
bsub -q $q -o $logdir/logCrossV_4.log CrossValidation 4  25 124 152 15 86 106           # Sample 4
bsub -q $q -o $logdir/logCrossV_5.log CrossValidation 5  25 93 113 15 100 122           # Sample 5
bsub -q $q -o $logdir/logCrossV_6.log CrossValidation 6  26 103 125 16 156 190          # Sample 6
bsub -q $q -o $logdir/logCrossV_7.log CrossValidation 7  19 78 95 11 143 174            # Sample 7
bsub -q $x -o $logdir/logCrossV_8.log CrossValidation 8  9 75 91 6 124 151      # Sample 8
bsub -q $q -o $logdir/logCrossV_9.log CrossValidation 9  26 92 112 16 140 171           # Sample 9
bsub -q $q -o $logdir/logCrossV_10.log CrossValidation 10  7 67 82 5 138 169            # Sample 10
bsub -q $q -o $logdir/logCrossV_11.log CrossValidation 11  25 121 148 15 130 159        # Sample 11
bsub -q $q -o $logdir/logCrossV_12.log CrossValidation 12  7 58 78 6 289 391 noRot      # Sample 12
bsub -q $q -o $logdir/logCrossV_13.log CrossValidation 13  7 53 71 6 221 299 noRot      # Sample 13
bsub -q $q -o $logdir/logCrossV_14.log CrossValidation 14  7 58 78 6 306 414 noRot      # Sample 14
bsub -q $q -o $logdir/logCrossV_15.log CrossValidation 15  7 66 90 6 149 201 noRot      # Sample 15
bsub -q $q -o $logdir/logCrossV_16.log CrossValidation 16  26 89 109 16 140 171         # Sample 16
bsub -q $q -o $logdir/logCrossV_17.log CrossValidation 17  26 75 92 16 190 233          # Sample 17
bsub -q $q -o $logdir/logCrossV_18.log CrossValidation 18  26 60 73 16 341 417          # Sample 18
bsub -q $q -o $logdir/logCrossV_19.log CrossValidation 19  26 54 66 16 275 336          # Sample 19
bsub -q $q -o $logdir/logCrossV_20.log CrossValidation 20  7 45 80 6 125 170 noRot      # Sample 20
bsub -q $q -o $logdir/logCrossV_21.log CrossValidation 21  26 73 90 16 266 325          # Sample 21
bsub -q $q -o $logdir/logCrossV_22.log CrossValidation 22  19 45 80 13 140 200 noRot    # Sample 22
bsub -q $q -o $logdir/logCrossV_23.log CrossValidation 23  26 72 89 16 148 180          # Sample 23
bsub -q $q -o $logdir/logCrossV_24.log CrossValidation 24  25 82 100 15 121 148         # Sample 24
bsub -q $q -o $logdir/logCrossV_25.log CrossValidation 25  26 62 76 16 220 269          # Sample 25
bsub -q $q -o $logdir/logCrossV_26.log CrossValidation 26  15 57 70 9 160 195           # Sample 26
bsub -q $q -o $logdir/logCrossV_27.log CrossValidation 27  26 63 76 16 182 222          # Sample 27
bsub -q $q -o $logdir/logCrossV_28.log CrossValidation 28  26 61 75 16 174 212          # Sample 28
bsub -q $q -o $logdir/logCrossV_29.log CrossValidation 29  26 56 68 16 179 219          # Sample 29

bsub -q $q -o $logdir/logCrossV_30.log CrossValidation 30  9 29 36 6 458 839            # Sample 30
bsub -q $q -o $logdir/logCrossV_31.log CrossValidation 31  25 29 35 15 470 574          # Sample 31
bsub -q $q -o $logdir/logCrossV_32.log CrossValidation 32  9 29 36 6 374 686            # Sample 32
bsub -q $q -o $logdir/logCrossV_33.log CrossValidation 33  25 21 25 15 740 904          # Sample 33
bsub -q $q -o $logdir/logCrossV_34.log CrossValidation 34  26 51 62 16 227 278          # Sample 34
bsub -q $q -o $logdir/logCrossV_35.log CrossValidation 35  25 52 64 15 321 392          # Sample 35
bsub -q $q -o $logdir/logCrossV_36.log CrossValidation 36  25 46 56 15 386 472          # Sample 36
bsub -q $q -o $logdir/logCrossV_37.log CrossValidation 37  25 54 66 15 287 350          # Sample 37
bsub -q $q -o $logdir/logCrossV_38.log CrossValidation 38  25 92 113 15 97 118          # Sample 38
bsub -q $q -o $logdir/logCrossV_39.log CrossValidation 39  25 88 108 15 93 114          # Sample 39
bsub -q $q -o $logdir/logCrossV_40.log CrossValidation 40  25 100 122 15 43 53          # Sample 40
bsub -q $q -o $logdir/logCrossV_41.log CrossValidation 41  25 69 84 15 96 117           # Sample 41
bsub -q $q -o $logdir/logCrossV_42.log CrossValidation 42  25 82 100 15 96 117          # Sample 42
bsub -q $q -o $logdir/logCrossV_43.log CrossValidation 43  25 125 153 15 50 61          # Sample 43
bsub -q $q -o $logdir/logCrossV_44.log CrossValidation 44  25 22 41 15 359 658          # Sample 44
bsub -q $q -o $logdir/logCrossV_45.log CrossValidation 45  26 70 86 16 155 285          # Sample 45
bsub -q $q -o $logdir/logCrossV_46.log CrossValidation 46  26 57 70 16 253 309          # Sample 46
bsub -q $q -o $logdir/logCrossV_47.log CrossValidation 47  25 100 122 15 130 159        # Sample 47
bsub -q $q -o $logdir/logCrossV_48.log CrossValidation 48  25 33 40 15 470 574          # Sample 48
bsub -q $q -o $logdir/logCrossV_49.log CrossValidation 49  25 19 24 15 525 642          # Sample 49
bsub -q $q -o $logdir/logCrossV_50.log CrossValidation 50  25 51 63 15 396 484          # Sample 50


