#!/bin/tcsh

cd ..
srtpath `cat .current` $BFARCH
cond24boot09
cd -

set q = xlong
set x = xxl
set d = `date +%y-%m-%d-%k`
set logdir = AWG82/logfiles/results/keys/CrossV/$d/
mkdir -p $logdir
rm $logdir/*

bsub -q $q -o $logdir/logCrossV_1.log CrossValidation 1         22 46 84        17 259 475      # Sample 1c
bsub -q $q -o $logdir/logCrossV_2.log CrossValidation 2         22 63 115       17 103 188      # Sample 2c
bsub -q $q -o $logdir/logCrossV_3.log CrossValidation 3         22 95 175       17 91 167       # Sample 3c
bsub -q $q -o $logdir/logCrossV_4.log CrossValidation 4         22 57 104       17 157 288      # Sample 4c
bsub -q $q -o $logdir/logCrossV_5.log CrossValidation 5         22 87 159       17 69 127       # Sample 5c
bsub -q $q -o $logdir/logCrossV_6.log CrossValidation 6         22 125 229      17 92 169       # Sample 6c
bsub -q $q -o $logdir/logCrossV_7.log CrossValidation 7         22 67 123       17 136 250      # Sample 7c
bsub -q $q -o $logdir/logCrossV_8.log CrossValidation 8         22 175 520      17 51 94        # Sample 8c
bsub -q $x -o $logdir/logCrossV_9.log CrossValidation 9         10 39 71        8 353 646       # Sample 9c
bsub -q $x -o $logdir/logCrossV_10.log CrossValidation 10       6 38 70         3 224 411       # Sample 10c
bsub -q $x -o $logdir/logCrossV_11.log CrossValidation 11       20 30 55        13 590 1082     # Sample 11c
bsub -q $x -o $logdir/logCrossV_12.log CrossValidation 12       10 44 80        8 150 275       # Sample 12c
bsub -q $x -o $logdir/logCrossV_13.log CrossValidation 13       10 45 83        8 111 204       # Sample 13c
bsub -q $q -o $logdir/logCrossV_14.log CrossValidation 14       22 32 59        17 359 1078     # Sample 14c
bsub -q $x -o $logdir/logCrossV_15.log CrossValidation 15       18 32 59        12 147 269      # Sample 15c
bsub -q $q -o $logdir/logCrossV_16.log CrossValidation 16       17 14 40        22 590 1082     # Sample 16c
bsub -q $q -o $logdir/logCrossV_17.log CrossValidation 17       22 45 83        17 136 249      # Sample 17c
bsub -q $q -o $logdir/logCrossV_18.log CrossValidation 18       22 32 59        17 196 360      # Sample 18c
bsub -q $q -o $logdir/logCrossV_19.log CrossValidation 19       22 48 88        17 144 264      # Sample 19c
bsub -q $q -o $logdir/logCrossV_20.log CrossValidation 20       22 38 70        17 229 419      # Sample 20c
bsub -q $x -o $logdir/logCrossV_21.log CrossValidation 21       22 20 37        14 570 1044     # Sample 21c
bsub -q $q -o $logdir/logCrossV_22.log CrossValidation 22       22 13 23        17 967 1773     # Sample 22c
bsub -q $q -o $logdir/logCrossV_23.log CrossValidation 23       22 18 34        17 713 1307     # Sample 23c
bsub -q $q -o $logdir/logCrossV_24.log CrossValidation 24       22 15 27        17 923 1693     # Sample 24c
bsub -q $q -o $logdir/logCrossV_25.log CrossValidation 25       22 32 59        17 286 524      # Sample 25c
bsub -q $q -o $logdir/logCrossV_26.log CrossValidation 26       22 43 79        17 239 439      # Sample 26c
bsub -q $q -o $logdir/logCrossV_27.log CrossValidation 27       22 28 51        17 290 532      # Sample 27c
bsub -q $q -o $logdir/logCrossV_28.log CrossValidation 28       22 54 98        17 127 233      # Sample 28c

bsub -q $q -o $logdir/logCrossV_29.log CrossValidation 29       22 32 59        17 286 524      # Sample 29c
bsub -q $q -o $logdir/logCrossV_30.log CrossValidation 30       22 43 79        17 239 439      # Sample 30c
bsub -q $q -o $logdir/logCrossV_31.log CrossValidation 31       22 28 51        17 290 532      # Sample 31c
bsub -q $q -o $logdir/logCrossV_32.log CrossValidation 32       22 54 98        17 127 233      # Sample 32c
bsub -q $q -o $logdir/logCrossV_41.log CrossValidation 41       22 34 63        17 195 357      # Sample 41c
bsub -q $q -o $logdir/logCrossV_42.log CrossValidation 42       22 39 71        17 158 290      # Sample 42c
bsub -q $q -o $logdir/logCrossV_43.log CrossValidation 43       22 64 117       17 103 188      # Sample 43c
bsub -q $q -o $logdir/logCrossV_44.log CrossValidation 44       22 47 86        17 158 290      # Sample 44c
bsub -q $q -o $logdir/logCrossV_45.log CrossValidation 45       22 72 133       12 103 188      # Sample 45c
bsub -q $q -o $logdir/logCrossV_46.log CrossValidation 46       13 111 204      7 51 93         # Sample 46c
bsub -q $q -o $logdir/logCrossV_47.log CrossValidation 47       22 82 150       17 93 171       # Sample 47c
bsub -q $q -o $logdir/logCrossV_48.log CrossValidation 48       22 80 147       17 66 121       # Sample 48c
bsub -q $q -o $logdir/logCrossV_51.log CrossValidation 51       22 58 107       17 124 227      # Sample 51c
bsub -q $q -o $logdir/logCrossV_52.log CrossValidation 52       22 95 175       17 38 70        # Sample 52c
bsub -q $q -o $logdir/logCrossV_53.log CrossValidation 53       22 39 72        17 170 312      # Sample 53c
bsub -q $q -o $logdir/logCrossV_54.log CrossValidation 54       22 139 254      17 93 170       # Sample 54c
bsub -q $q -o $logdir/logCrossV_55.log CrossValidation 55       22 69 127       17 117 214      # Sample 55c
bsub -q $q -o $logdir/logCrossV_56.log CrossValidation 56       22 126 232      17 45 82        # Sample 56c
bsub -q $q -o $logdir/logCrossV_57.log CrossValidation 57       22 49 90        17 130 238      # Sample 57c
bsub -q $q -o $logdir/logCrossV_58.log CrossValidation 58       22 155 284      17 64 117       # Sample 58c
bsub -q $q -o $logdir/logCrossV_61.log CrossValidation 61       22 58 107       17 124 227      # Sample 61c
bsub -q $q -o $logdir/logCrossV_62.log CrossValidation 62       22 95 175       17 38 70        # Sample 62c
bsub -q $q -o $logdir/logCrossV_63.log CrossValidation 63       22 39 72        17 170 312      # Sample 63c
bsub -q $q -o $logdir/logCrossV_64.log CrossValidation 64       22 139 254      17 93 170       # Sample 64c
bsub -q $q -o $logdir/logCrossV_65.log CrossValidation 65       22 69 127       17 117 214      # Sample 65c
bsub -q $q -o $logdir/logCrossV_66.log CrossValidation 66       22 126 232      17 45 82        # Sample 66c
bsub -q $q -o $logdir/logCrossV_67.log CrossValidation 67       22 49 90        17 130 238      # Sample 67c
bsub -q $q -o $logdir/logCrossV_68.log CrossValidation 68       22 155 284      17 64 117       # Sample 68c

